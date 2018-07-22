#!/usr/bin/env python

"""
NAME
    elas.py

BACKGROUND
    This is a python implementation of a fortran program by Rebecca Brannon
    (brannon@mech.utah.edu)

PURPOSE
    Given any two INDEPENDENT elastic constants (from the list below) this
    program will compute all of the others. This program also contains formulas
    that will give you an estimate of how the elastic properties will change when
    the material contains pores or cracks. If you don't care about this
    information, simply enter zero (0) when asked to provide a damage measure (or
    just hit return because undamaged properties are the default). If you ARE
    interested in getting the changes in elastic properties that result from
    PORES, then first enter any two independent elastic properties of the
    nonporous matrix material, and then enter the pore volume fraction when asked
    for the damage measure. If you seek the changes in elastic properties
    resulting from CRACKS, then specify the NEGATIVE of the number of cracks per
    unit volume times the mean of the cube of the crack radius (using a negative
    is simply a low-tech way of telling this program that the damage measure
    represents crack data instead of porosity data). Again note: the first two
    parameters must be the moduli of the undamaged matrix material. Solving the
    inverse problem of getting the undamaged moduli when the damaged moduli are
    known requires an iterative solver -- contact me (brannon@mech.utah.edu) if
    you need routines that will do that.

USAGE
    Interactive

     % elas

    elas will ask to enter values of elastic constants, as follows

       FIRST elastic constant? E = 200e9
       SECOND elastic constant? G = 79e9
       Density? [None] (enter zero if wave speeds not desired)

   Defaults [in brackets] may be accepted by hitting return. After entering data
   one time, the program will loop back for a new set of inputs, and all your old
   inputs will be the new defaults (very useful for correcting typos).

   Direct

     % elas --emod1=val --emod2=val

   Execute elas -h for possible --emod[1,2] strings

   A NOTE ABOUT WAVE SPEEDS...
   When shock physicists say `cs', they mean what we call `cl'. Longitudinal
   waves (called "P-waves" by geologists) travel at speed `cl'. Shear waves
   ("S-waves") travel at speed ct. When you compute elastic moduli using wave
   speeds, be aware that you are actually finding the ISENTROPIC moduli.
   Conversely, if you seek the wave speeds, then you must enter ISENTROPIC
   moduli. For isotropic materials, converting between isentropic moduli and
   isothermal moduli is straightforward, but this program won't do it for you.
   Wave speeds correspond to ISENTROPIC moduli because sound waves travel such
   that the stress-strain response is the same as that when the entropy is
   constant. elastic moduli that are measured quasistatically using a load cell
   are typically determined from stress-strain response at constant temperature
   (e.g., room temperature) Here is how the conversions are computed... Let
   subscripts "s" and "T" indicate isentropic and isothermal, respectively. Then
   relevant formulas are: G_s = G_T and K_s/K_T = c_p/c_v, where c_p and c_v are
   specific heats at const. pressure and volume. If you don't have values for
   both types of specific heats, consider looking at
   http://mech.utah.edu/~brannon/DerivativeRecursionTables.pdf for alternative
   formulas.
"""

constformulas="""
   #########################################################################
   Lame = First Lame parameter         = G(E-2G)/(3G-E) = K-2G/3
      G = Shear modulus (= 2nd Lame parameter, mu)   = E/2/(1+Nu)
      E = Young's modulus              = 3K(1-2Nu) = 2G(1+Nu) = 9KG/(3K+G)
     Nu = Poisson's ratio              = (3K-E)/6K = Lame/2/(Lame+G)
      K = bulk modulus                 = E/3/(1-2Nu) = Lame + 2G/3
      H = constrained modulus          = 2G+Lame = 3K-2Lame = K + 4G/3
     Ko = SIGy/SIGx in uniaxial strain = Nu/(1-Nu)
     Cl = longitudinal wave speed      = Sqrt[H/Rho]
     Ct = shear (TRANSVERSE) wave speed  = Sqrt[G/Rho)
     Co = bulk/plastic wave speed      = Sqrt[K/Rho]=Sqrt[Cl^2-4(Ct^2)/3]
     Cr = thin rod elastic wave speed  = Sqrt[E/Rho]
   #########################################################################
"""
import sys
import os
import re
import argparse
import readline
from math import sqrt

LAME, G, E, NU, K, H = 'Lame', 'G', 'E', 'Nu', 'K', 'H'
C10, D1 = 'C10', 'D1'
KO, CL, CT, CO, CR = 'Ko', 'Cl', 'Ct', 'Co', 'Cr'
RHO = 'Rho'

exe = 'elas'
manpage = '{0} \n {1}'.format(__doc__, constformulas)
interactive_help = ('Enter elastic constants in the form\n'
                    '\tconst = val\n\n{0} recognizes \n {1}\n'
                    'Enter h to display this message again\n'
                    'Enter q to quit\n'.format(exe, constformulas))

class ElasticConstantsError(Exception):
    pass

class NonPositiveDefiniteError(Exception):
    def __init__(self, postscript='dum'):
        message = 'nonpositive definite elastic constants,' + postscript
        super(NonPositiveDefiniteError, self).__init__(message)

class AmbiguousInputError(Exception):
    def __init__(self):
        message = 'Ambiguous elastic constants'
        super(AmbiguousInputError, self).__init__(message)

class ElasticConstants(object):
    def __init__(self):
        self.names = (LAME, G, E, NU, K, H, KO, CL, CT, CO, CR, RHO)
        self.nconsts = len(self.names)
        self.constants = dict(zip(self.names, [None]*self.nconsts))
        self.moduli = self.names[:7]
        self.wavespeeds = self.names[7:-1]

    def __contains__(self, item):
        return item.capitalize() in self.names

    def convert_aliases(self, **kwargs):
        kwds = {}
        for (kwd, val) in kwargs.items():
            if kwd.lower() == 'mu':
                kwds[G] = val
            else:
                kwds[kwd] = val
        return kwds

    def compute_elastic_constants(self, **kwargs):
        """Given any two elastic elastic, compute all remaining constants

        """
        kwargs = self.convert_aliases(**kwargs)
        kwds = dict([(key.capitalize(), value) for (key, value) in kwargs.items()])
        calledwith = ' Called with: ' + repr(kwds)

        # check for goodness of input
        rho = kwds.pop(RHO, None)
        if rho is not None and rho <= 0.:
            raise ElasticConstantsError('Expected RHO to be > 0, got '
                                        '{0}'.format(rho))
        if kwds.get(NU) is not None:
            if kwds[NU] < -1 or kwds[NU] >= .5:
                raise ElasticConstantsError('Expected -1 < NU <= .5, got '
                                            '{0}'.format(kwds[NU]))
        if len(kwds) != 2:
            raise ElasticConstantsError('Expected 2 elastic constants, '
                                        'got {0}'.format(len(kwds)))

        # determine indices of passed constants
        ij = []
        errors = []
        for (key, value) in kwds.items():
            try:
                ij.append(self.names.index(key))
            except ValueError:
                errors.append(key)
        if errors:
            raise ElasticConstantsError('Unexpected elastic constants: '
                                        '{0}'.format(', '.join(errors)))
        ii, jj = ij

        # determine if density is required
        kwds[RHO] = rho
        if (ii > 6 and jj > 6) and not kwds[RHO]:
            raise ElasticConstantsError('Density must be given when a wave '
                                        'speed is given')

        if ii == 6: # ko
            ii = 3
            kwds[NU] = kwds[KO] / (1. + kwds[KO])
        elif ii == 7: # cl
            ii = 5
            kwds[H] = kwds[RHO] * kwds[CL] * kwds[CL]
        elif ii == 8: # ct
            ii = 1
            kwds[G] = kwds[RHO] * kwds[CT] * kwds[CT]
        elif ii == 9: # co
            ii = 4
            kwds[K] = kwds[RHO] * kwds[CO] * kwds[CO]
        elif ii == 10: # cr
            ii = 2
            kwds[E] = kwds[RHO] * kwds[CR] * kwds[CR]

        if jj == 6: # ko
            jj = 3
            kwds[NU] = kwds[KO] / (1. + kwds[KO])
        elif jj == 7: # cl
            jj = 5
            kwds[H] = kwds[RHO] * kwds[CL] * kwds[CL]
        elif jj == 8: # ct
            jj = 1
            kwds[G] = kwds[RHO] * kwds[CT] * kwds[CT]
        elif jj == 9: # co
            jj = 4
            kwds[K] = kwds[RHO] * kwds[CO] * kwds[CO]
        elif jj == 10: # cr
            jj = 2
            kwds[E] = kwds[RHO] * kwds[CR] * kwds[CR]

        #  At this point, ii and jj each range from 1 to 6, and are distinct.
        #  There are 15 possible ways to choose 2 numbers from 6:
        case = 0
        if ii < jj:
            case = 10 * (ii + 1) + jj + 1
        else:
            case = 10 * (jj + 1) + ii + 1

        # For each case, determine G and NU. From these two, all other
        # constants will be determined below
        if case == 12: # Lame, G
            if kwds[LAME] + kwds[G] == 0.:
                raise NonPositiveDefiniteError('LAME+G'+calledwith)
            kwds[NU] = kwds[LAME] / 2. / (kwds[LAME] + kwds[G])
            if abs(kwds[G]) < abs(kwds[LAME])*1.e-9:
                kwds[G] = max(0.0, kwds[G])
        elif case == 13: # Lame, E
            A = kwds[E] * kwds[E]
            A += 2. * kwds[LAME] * kwds[E]
            A += 9. * kwds[LAME] * kwds[LAME]
            if A < 0.:
                raise NonPositiveDefiniteError('LAME+E, A<0'+calledwith)
            A = sqrt(A)
            kwds[G] = (A - 3. * kwds[LAME] + kwds[E]) / 4.
            kwds[NU] = (A - kwds[E] - kwds[LAME]) / 4. / kwds[LAME]
        elif case == 14: # Lame, nu
            if kwds[NU] == 0.:
                raise AmbiguousInputError
            kwds[G]  = kwds[LAME] * (1. - 2. * kwds[NU]) / 2. / kwds[NU]
        elif case == 15: # Lame,K
            if 3. * kwds[K] - kwds[LAME] == 0.:
                raise NonPositiveDefiniteError('LAME+K'+calledwith)
            kwds[G]  = 3. * (kwds[K] - kwds[LAME]) / 2.
            kwds[NU] = kwds[LAME] / (3. * kwds[K] - kwds[LAME])
        elif case == 16: # Lame, H
            if kwds[H] + kwds[LAME] == 0.:
                raise NonPositiveDefiniteError('LAME+H'+calledwith)
            kwds[G]  = (kwds[H] - kwds[LAME]) / 2.
            kwds[NU] = kwds[LAME] / (kwds[H] + kwds[LAME])
        elif case == 23: # G, E
            kwds[NU] = (kwds[E] - 2. * kwds[G]) / 2. / kwds[G]
        elif case == 24: # G, nu
            pass
        elif case == 25: # G,K
            if 3. * kwds[K] + kwds[G] == 0.:
                raise NonPositiveDefiniteError('G+K'+calledwith)
            kwds[NU] = (3.*kwds[K]-2.*kwds[G])/2./(3.*kwds[K]+kwds[G])
        elif case == 26: # G, H
            if kwds[H] - kwds[G] == 0.:
                raise NonPositiveDefiniteError('G+H'+calledwith)
            kwds[NU] = (kwds[H]-2.*kwds[G])/2./(kwds[H]-kwds[G])
        elif case == 34: # E, nu
            if 1. + kwds[NU] == 0.:
                raise NonPositiveDefiniteError('E+NU'+calledwith)
            kwds[G]  = kwds[E] / 2. / (1. + kwds[NU])
        elif case == 35: # E,K
            if 9. * kwds[K] - kwds[E] == 0.  or  kwds[K] == 0.:
                raise NonPositiveDefiniteError('E+K'+calledwith)
            kwds[G]  = 3. * kwds[E] * kwds[K] / (9. * kwds[K] - kwds[E])
            kwds[NU] = (3. * kwds[K] - kwds[E]) / 6. / kwds[K]
        elif case == 36: # E, H
            B = kwds[E] * kwds[E]
            B += 9. * kwds[H] * kwds[H]
            B -= 10. * kwds[E] * kwds[H]
            if B <0.  or  kwds[H] == 0.:
                raise NonPositiveDefiniteError('E+H'+calledwith)
            B = sqrt(B)
            kwds[G]  = (3. * kwds[H] - B + kwds[E]) / 8.
            kwds[NU] = (B - kwds[H] + kwds[E]) / 4. / kwds[H]
        elif case == 45: # nu,K
            if 1. + kwds[NU] == 0.:
                raise NonPositiveDefiniteError('NU+K'+calledwith)
            kwds[G]  = 3.*kwds[K]*(1.-2.*kwds[NU])/2./(1.+kwds[NU])
        elif case == 46: # nu, H
            if 1. - kwds[NU] == 0.:
                raise NonPositiveDefiniteError('NU+H'+calledwith)
            kwds[G]  = kwds[H]*(1.-2.*kwds[NU])/2./(1.-kwds[NU])
        elif case == 56: # K, H
            if 3. * kwds[K] + kwds[H] == 0.:
                raise NonPositiveDefiniteError('K+H'+calledwith)
            kwds[G]  = 3. * (kwds[H] - kwds[K]) / 4.
            kwds[NU] = (3.*kwds[K]-kwds[H])/(3.*kwds[K]+kwds[H])
        else:
            raise ElasticConstantsError('Unexpected case')

        kwds[LAME] = 2. * kwds[G] * kwds[NU] / (1. - 2. * kwds[NU])
        kwds[E] = 2. * kwds[G] * (1. + kwds[NU])
        kwds[K] = 2.*kwds[G] * (1. + kwds[NU]) / 3. / (1. - 2. * kwds[NU])
        kwds[H] = 2. * kwds[G] * (1. - kwds[NU]) / (1 - 2 * kwds[NU])
        kwds[KO] = kwds[NU] / (1. - kwds[NU])

        if kwds[G] < 0. or kwds[K] < 0.:
            message = ('K and G must both be non-negative: '
                       'K={0:.4e} G={1:.4e}'.format(kwds[K],kwds[G])+calledwith)
            raise NonPositiveDefiniteError(message)

        if kwds[RHO]:
            kwds[CL] = sqrt(kwds[H] / kwds[RHO])
            kwds[CT] = sqrt(kwds[G] / kwds[RHO])
            kwds[CO] = sqrt(kwds[K] / kwds[RHO])
            kwds[CR] = sqrt(kwds[E] / kwds[RHO])

        kwds[C10] = kwds[E] / (4. * (1. + kwds[NU]))
        kwds[D1] = 6. * (1. - 2. * kwds[NU]) / kwds[E]

        self.constants.update(kwds)
        return kwds

    def print_elastic_constants(self):
        if not self.constants:
            raise ElasticConstantsError('constants not yet computed')

        print('\nElastic moduli:')
        for name in self.moduli:
            print('{0} = {1:g}'.format(name, self.constants[name]))
            continue

        print('\nNeo-Hooke Hyperelastic constants:')
        print('{0} = {1:g}'.format(C10, self.constants[C10]))
        print('{0} = {1:g}'.format(D1, self.constants[D1]))

        if self.constants.get(RHO):
            print('\nWavespeeds:')
            print('Density = {0:g}'.format(self.constants[RHO]))
            for name in self.wavespeeds:
                print('{0} = {1:g}'.format(name, self.constants[name]))
                continue
        print('')
        return

    def interactive_guide(self):
        """Command line interface"""
        print(interactive_help)
        aa = [None, None]
        bb = [None, None]
        rr = [RHO, None]

        while True:
            self.get_input('First elastic constant? [default {0}={1}]', aa)
            while True:
                self.get_input('Second elastic constant? [default {0}={1}]', bb)

                if aa[0] == bb[0]:
                    print('Second elastic constant must differ '
                          'from first\nTry again')
                    continue
                break

            self.get_input('Density? [default {1}]', rr, optional=1)

            ui = {}
            ui[aa[0]] = aa[1]
            ui[bb[0]] = bb[1]
            ui[rr[0]] = rr[1]
            try:
                self.compute_elastic_constants(**ui)
                self.print_elastic_constants()
            except (NonPositiveDefiniteError, AmbiguousInputError) as e:
                message = '{0}, try again\n'.format(e.args[0])
                print(message)

            continue

        return

    def get_input(self, query, args, optional=0):

        while True:

            # ask for name, value pairs
            try:
                o = str(input('{0}: '.format(query.format(*args))))
            except (KeyboardInterrupt, EOFError):
                raise SystemExit('\n')

            if not o.split() and optional:
                return

            inp = [x.strip().capitalize()
                   for x in re.split(r'[=,]', o) if x.split()]

            if not inp:
                inp = [x for x in args]

            elif inp[0] == 'Q':
                raise SystemExit('done')

            elif inp[0] == 'H':
                print(interactive_help)
                continue

            if len(inp) == 1:
                # use default name
                inp.insert(0, args[0])

            if optional and inp[1] is None:
                return

            try:
                name, val = inp[0], float(inp[1])
            except ValueError:
                print('\nExpected name = value, got {0}\nTry again\n'.format(o))
                continue

            if name not in self.names:
                print('\nUnknown constant {0}\nvalid entries are '
                      '{1}\nTry again'.format(name, ', '.join(self.names)))
                continue

            break

        args[:2] = [name, val]

        return


def main(argv=None):
    """ Fetches user input and sends in right format to compute_elastic_constants

    """
    argv = argv or sys.argv[1:]

    # -- command line option parsing
    p= argparse.ArgumentParser()
    p.add_argument('--lam', dest=LAME, help='First Lame parameter')
    p.add_argument('--shear', dest=G, help='Shear modulus')
    p.add_argument('--youngs', dest=E, help="Young's modulus")
    p.add_argument('--poissons', dest=NU, help="Poisson's ratio")
    p.add_argument('--bulk', dest=K, help='Bulk modulus')
    p.add_argument('--constrained', dest=H, help='Constrained modulus')
    p.add_argument('--ko', dest=KO, help='SIGy/SIGx in uniaxial strain')
    p.add_argument('--cl', dest=CL, help='Longitudinal wave speed')
    p.add_argument('--ct', dest=CT, help='Shear (TRANSVERSE) wave speed')
    p.add_argument('--co', dest=CO, help='Bulk/plastic wave speed')
    p.add_argument('--cr', dest=CR, help='Thin rod elastic wave speed')
    p.add_argument('--rho', dest=RHO, help='Density')
    p.add_argument('--man', action='store_true', help='display manpage')
    args = p.parse_args(argv)

    if args.man:
        parser.print_help()
        sys.exit(manpage)

    EC = ElasticConstants()

    # check if user gave input
    ui = {}
    for name in EC.names:
        try: ui[name] = float(getattr(args, name))
        except (AttributeError, TypeError): continue

    if ui:
        EC.compute_elastic_constants(**ui)
        EC.print_elastic_constants()
        return 0

    EC.interactive_guide()
    return 0

def elas(**kwargs):
    EC = ElasticConstants()
    return EC.compute_elastic_constants(**kwargs)

if __name__ == '__main__':
    main()
