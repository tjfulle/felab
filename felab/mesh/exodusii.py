import re
import os
import sys
import logging
import datetime
from numpy import *
from os.path import basename, join, splitext, isfile

from ._mesh import element_block
from ..constants import *
from .post import StepRepository1
from ..elemlib import element_family

# True if we are running on Python 3.
PY3 = sys.version_info[0] == 3
try:
    NETCDF4 = True
    from netCDF4 import Dataset
except ImportError:
    NETCDF4 = False
    warn_all = False
    if PY3 and warn_all:
        logging.warn('USING SCIPY.IO.NETCDF LIBRARY TO WRITE EXODUSII FILES.  '
                     'INSTALL NETCDF4 FOR BETTER RELIABILITY.')
    from .netcdf import netcdf_file as Dataset

__all__ = ['File', 'put_nodal_solution']

def cat(*args):
    return ''.join('{0}'.format(a).strip() for a in args)

def stringify2(a):
    try:
        a = a.tolist()
    except AttributeError:
        pass
    try:
        return [''.join(b.decode() for b in row if b.split()) for row in a]
    except AttributeError:
        return ['']

def stringify(a):
    return ''.join(a.decode().strip())

def aindex(arr, ij):
    a = arr.tolist()
    return array([a.index(i) for i in ij])

def sortexoname(key):
    if len(key) == 1:
        return (ord(key), 20)
    s = key.lower()
    if s.startswith('displ'):
        return (0, 'xyz'.index(s[-1]))
    if s.endswith('xx'):
        return (key[:-2], 0)
    elif s.endswith('yy'):
        return (key[:-2], 1)
    elif s.endswith('zz'):
        return (key[:-2], 3)
    elif s.endswith('xy'):
        return (key[:-2], 4)
    elif s.endswith('yz'):
        return (key[:-2], 5)
    elif s.endswith('xz'):
        return (key[:-2], 6)
    elif s.endswith('x'):
        return ord(key[0]), 0
    elif s.endswith('y'):
        return ord(key[0]), 1
    elif s.endswith('z'):
        return ord(key[0]), 2
    else:
        return ord(key[0]), 10

def labels_and_data(key, data, numdim, numnod=None):
    key = key.upper()
    if key == 'U':
        # displacement is a special case
        labels = VAR_DISP_NAMES(numdim)
        if data.ndim == 1:
            data = data.reshape(-1, 1)
        return labels, data
    if data.ndim == 1:
        labels = [key]
        data = data.reshape(-1, 1)
    elif key in ('T',):
        # Scalars
        assert data.shape[1] == 1
        labels = [key]
    elif key in ('R',):
        # Vectors
        labels = [key + x for x in 'xyz'[:numdim]]
    elif key in ('S', 'E', 'DE'):
        ndir, nshr = {6: (3,3), 4: (3, 1), 3: (2, 1)}[data.shape[-1]]
        # Symmetric tensors
        components = ('xx','yy','zz')[:ndir] + ('xy','yz','xz')[:nshr]
        labels = [key + x for x in components]
    else:
        raise KeyError(key)
    return tuple(labels), data

def stringtoarr(string, NUMCHARS):
    # function to convert a string to a array of NUMCHARS characters
    arr = zeros(NUMCHARS, 'S1')
    arr[0:len(string)] = tuple(string)
    return arr

def adjstr(string):
    return '{0:{1}s}'.format(string, LEN_STRING)[:LEN_STRING]

# --- Data types
FLOAT = 'f'
CHAR = 'c'
INT = 'i'

# --- Global dimensions and variables
DIM_LEN_STRING     = 'len_string'
LEN_STRING         = 33
DIM_LEN_LINE       = 'len_line'
DIM_FOUR           = 'four'
DIM_NUM_DIM        = 'num_dim'
DIM_NUM_QA         = 'num_qa_rec'
DIM_TIME_STEP      = 'time_step'
DIM_MAX_STEPS      = 'max_steps'
VAR_TIME_WHOLE     = 'time_whole'
VAR_QA_RECORDS     = 'qa_records'
VAR_COOR_NAMES     = 'coor_names'
VAR_COOR_NAME      = lambda i: cat('coord', 'xyz'[i])
VAR_DISP_NAMES     = lambda n: [cat('displ', 'xyz'[i]) for i in range(n)]

# --- Element dimensions and variables
DIM_NUM_ELE        = 'num_elem'
DIM_NUM_ELE_VAR    = 'num_elem_var'
DIM_NUM_ELEBLK     = 'num_el_blk'
VAR_NAME_ELE_VAR   = 'name_elem_var'
VAR_EB_PROP1       = 'eb_prop1'
VAR_EB_STATUS      = 'eb_status'
VAR_EB_NAMES       = 'eb_names'
VAR_ELE_MAP        = lambda i: cat('elem_map', i)
VAR_ELE_NUM_MAP    = 'elem_num_map'
VAR_ELE_TAB        = 'elem_var_tab'
DIM_NUM_EL_IN_BLK  = lambda i: cat('num_el_in_blk', i)
DIM_NUM_NOD_PER_EL = lambda i: cat('num_nod_per_el', i)
VAR_BLKCON         = lambda i: cat('connect', i)
VALS_ELE_VAR       = lambda i, j: cat('vals_elem_var', i, 'eb', j)
DIM_NUM_ES         = 'num_elem_sets'
DIM_NUM_ELE_ES     = lambda i: cat('num_ele_els', i)
VAR_ES_PROP1       = 'els_prop1'
VAR_ES_NAMES       = 'els_names'
VAR_ELE_ES         = lambda i: cat('elem_els', i)


# --- Node dimensions and variables
DIM_NUM_NOD        = 'num_nodes'
DIM_NUM_NOD_VAR    = 'num_nod_var'
VAR_NOD_MAP        = lambda i: cat('node_map', i)
VAR_NOD_NUM_MAP    = 'node_num_map'
VAR_NAME_NOD_VAR   = 'name_nod_var'
VALS_NOD_VAR       = lambda i: cat('vals_nod_var', i)

# --- Node set dimensons and variables
DIM_NUM_NS         = 'num_node_sets'
DIM_NUM_NOD_NS     = lambda i: cat('num_nod_ns', i)
VAR_NS_PROP1       = 'ns_prop1'
VAR_NS_NAMES       = 'ns_names'
VAR_NOD_NS         = lambda i: cat('node_ns', i)

# --- Global variable dimensions and variables
DIM_NUM_GLO_VAR    = 'num_glo_var'
VALS_GLO_VAR       = 'vals_glo_var'

# --- Side set dimensions and variables
DIM_NUM_SS         = 'num_side_sets'
DIM_NUM_SIDE_SS    = lambda i: cat('num_side_ss', i)
DIM_NUM_ELE_SS     = lambda i: cat('num_elem_ss', i)
VAR_SS_PROP1       = 'ss_prop1'
VAR_SS_NAMES       = 'ss_names'
VAR_SIDE_SS        = lambda i: cat('side_ss', i)
VAR_ELE_SS         = lambda i: cat('elem_ss', i)

# --- Field variable dimensions and variables (femlib specific)
DIM_NUM_FIELD      = 'num_field'
VAR_STEP_NUM       = 'step_num'
VAR_STEP_NAMES     = 'step_names'
VAR_FIELD_ELE_VAR  = 'field_elem_var'
VAR_FIELD_NOD_VAR  = 'field_nod_var'
VAR_FO_PROP1       = 'fo_prop1'
VAR_FO_NAMES       = 'fo_names'
VAR_FO_TYPES       = 'fo_types'
VAR_FO_VALINV      = 'fo_valinv'

def File(filename, mode='r'):
    if mode not in 'rw':
        raise ValueError('UNKNOWN FILE MODE {0}'.format(mode))
    if mode == 'r':
        return EXOFileReader(filename)
    return EXOFileWriter(filename)

def get_exo_eletyp(numdim, numnod):
    if numnod == 2:
        eletyp = 'TRUSS'
    elif numdim == 1:
        eletyp = 'TRUSS'
    elif numdim == 2:
        if numnod in (3, 6):
            eletyp = 'TRI'
        elif numnod in (4, 8):
            eletyp = 'QUAD'
    elif numdim == 3:
        if numnod in (4, 6):
            eletyp = 'TET'
        elif numnod in (8, 20):
            eletyp = 'HEX'
    else:
        raise ValueError('UNKNOWN ELEMENT TYPE')
    return eletyp

class EXOFile(object):
    mode = None
    def close(self):
        self.fh.close()

    def getdim(self, name, default=None):
        try:
            x = self.fh.dimensions[name]
        except KeyError:
            return default
        if NETCDF4:
            dim = x.size
        else:
            dim = x
        if dim is None and name != DIM_TIME_STEP:
            return 0
        return int(dim)

    def setncattr(self, variable, name, value):
        if NETCDF4:
            self.fh.variables[variable].setncattr(name, value)
        else:
            setattr(self.fh.variables[variable], name, value)

    def open_file(self, filename, mode='r'):
        if NETCDF4:
            fh = Dataset(filename, mode=mode, format='NETCDF4_CLASSIC')
        else:
            fh = Dataset(filename, mode=mode)
        return fh

class EXOFileWriter(EXOFile):
    mode = 'w'
    def __init__(self, filename):
        '''
        Notes
        -----
        The EXOFile class is an interface to the Exodus II api. Its methods
        are named after the analogous method from the Exodus II C bindings,
        minus the prefix 'ex_'.

        '''
        self.fh = self.open_file(filename, mode='w')
        self.jobid = splitext(basename(filename))[0]
        self.filename = filename
        self.initialized = False
        self.count = 0

    def update(self):
        pass

    def genesis(self, nodmap, elemap, coord, element_blocks,
                nodesets=None, elemsets=None, sidesets=None):

        numele = len(elemap)
        numnod = coord.shape[0]
        if coord.ndim == 1:
            coord = coord.reshape(-1, 1)
        numnod, numdim = coord.shape
        # Node map: nodmap1[i] is the external node label of the ith node
        nodmap1 = array(sorted(nodmap.keys(), key=lambda k: nodmap[k]))

        elemap = elemap
        elemap1 = array(sorted(elemap.keys(), key=lambda k: elemap[k]))

        # ------------------------------------------------------------------- #
        # ------------------------------------------------- record arrays --- #
        # ------------------------------------------------------------------- #
        self.fh.createDimension(DIM_TIME_STEP, None)
        self.fh.createVariable(VAR_TIME_WHOLE, FLOAT, (DIM_TIME_STEP,))

        # ------------------------------------------------------------------- #
        # -------------------------------- standard ExodusII dimensioning --- #
        # ------------------------------------------------------------------- #
        self.fh.floating_point_word_size = 4
        self.fh.version = 5.0300002
        self.fh.file_size = 1
        self.fh.api_version = 5.0300002
        self.fh.title = 'finite element simulation'

        self.fh.filename = basename(self.filename)
        self.fh.jobid = self.jobid

        self.fh.createDimension(DIM_LEN_STRING, LEN_STRING)
        self.fh.createDimension(DIM_LEN_LINE, 81)
        self.fh.createDimension(DIM_FOUR, 4)

        self.fh.createDimension(DIM_NUM_DIM, numdim)
        self.fh.createDimension(DIM_NUM_NOD, numnod)
        self.fh.createDimension(DIM_NUM_ELE, numele)

        # node and element number maps
        self.fh.createVariable(VAR_NOD_MAP(1), INT, (DIM_NUM_NOD,))
        self.fh.variables[VAR_NOD_MAP(1)][:] = nodmap1

        # ------------------------------------------------------------------- #
        # ---------------------------------------------------- QA records --- #
        # ------------------------------------------------------------------- #
        now = datetime.datetime.now()
        day = now.strftime("%m/%d/%y")
        hour = now.strftime("%H:%M:%S")
        self.fh.createDimension(DIM_NUM_QA, 1)
        self.fh.createVariable(VAR_QA_RECORDS, CHAR,
                               (DIM_NUM_QA, DIM_FOUR, DIM_LEN_STRING))
        self.fh.variables[VAR_QA_RECORDS][0, 0, :] = adjstr('felab')
        self.fh.variables[VAR_QA_RECORDS][0, 1, :] = adjstr(self.jobid)
        self.fh.variables[VAR_QA_RECORDS][0, 2, :] = adjstr(day)
        self.fh.variables[VAR_QA_RECORDS][0, 3, :] = adjstr(hour)

        # ------------------------------------------------------------------- #
        # ------------------------------------------ node coordinate data --- #
        # ------------------------------------------------------------------- #
        self.fh.createVariable(VAR_COOR_NAMES, CHAR, (DIM_NUM_DIM, DIM_LEN_STRING))
        for i in range(numdim):
            self.fh.variables[VAR_COOR_NAMES][i,:] = adjstr(VAR_COOR_NAME(i))
            self.fh.createVariable(VAR_COOR_NAME(i), FLOAT, (DIM_NUM_NOD,))
            self.fh.variables[VAR_COOR_NAME(i)][:] = coord[:, i]

        # ------------------------------------------------------------------- #
        # --------------------------------------- element block meta data --- #
        # ------------------------------------------------------------------- #
        # block IDs - standard map
        # element map
        num_el_blk = len(element_blocks)
        self.fh.createDimension(DIM_NUM_ELEBLK, num_el_blk)

        self.fh.createVariable(VAR_EB_PROP1, INT, (DIM_NUM_ELEBLK,))
        self.setncattr(VAR_EB_PROP1, 'name', 'ID')

        self.fh.createVariable(VAR_ELE_MAP(1), INT, (DIM_NUM_ELE,))
        self.fh.variables[VAR_ELE_MAP(1)][:] = elemap1

        self.fh.createVariable(VAR_EB_STATUS, INT, (DIM_NUM_ELEBLK,))
        self.fh.variables[VAR_EB_STATUS][:] = ones(num_el_blk, dtype=int)

        self.fh.createVariable(VAR_EB_NAMES, CHAR, (DIM_NUM_ELEBLK, DIM_LEN_STRING))
        for (ieb, eb) in enumerate(element_blocks, start=1):
            self.fh.variables[VAR_EB_NAMES][ieb-1,:] = adjstr(eb.name)
            self.fh.variables[VAR_EB_PROP1][ieb-1] = eb.id

            # block connect
            blkcon = eb.elecon + 1

            ne, nn = blkcon.shape
            d1, d2 = DIM_NUM_EL_IN_BLK(ieb), DIM_NUM_NOD_PER_EL(ieb)
            self.fh.createDimension(d1, ne)
            self.fh.createDimension(d2, nn)

            # set up the element block connectivity
            self.fh.createVariable(VAR_BLKCON(ieb), INT, (d1, d2))
            self.fh.variables[VAR_BLKCON(ieb)][:] = blkcon

            # element type
            elem_type = get_exo_eletyp(numdim, nn)
            self.fh.variables[VAR_BLKCON(ieb)].elem_type = elem_type

        # ------------------------------------------------------------------- #
        # ----------------------------------------- element set meta data --- #
        # ------------------------------------------------------------------- #
        if elemsets:
            nes = len(elemsets)
            self.fh.createDimension(DIM_NUM_ES, nes)
            prop1 = arange(nes, dtype=int32) + 1
            self.fh.createVariable(VAR_ES_PROP1, INT, (DIM_NUM_ES,))
            self.fh.variables[VAR_ES_PROP1][:] = prop1
            self.setncattr(VAR_ES_PROP1, 'name', 'ID')
            self.fh.createVariable(VAR_ES_NAMES, CHAR,
                                   (DIM_NUM_ES, DIM_LEN_STRING))
            i = 1
            for (name, elems) in elemsets.items():
                self.fh.variables[VAR_ES_NAMES][i-1,:] = adjstr(name)
                self.fh.createDimension(DIM_NUM_ELE_ES(i), len(elems))
                self.fh.createVariable(VAR_ELE_ES(i), INT, (DIM_NUM_ELE_ES(i),))
                self.fh.variables[VAR_ELE_ES(i)][:] = elems + 1
                i += 1

        # ------------------------------------------------------------------- #
        # -------------------------------------------- node set meta data --- #
        # ------------------------------------------------------------------- #
        if nodesets:
            nns = len(nodesets)
            self.fh.createDimension(DIM_NUM_NS, nns)
            # node set IDs - standard map
            prop1 = arange(nns, dtype=int32) + 1
            self.fh.createVariable(VAR_NS_PROP1, INT, (DIM_NUM_NS,))
            self.fh.variables[VAR_NS_PROP1][:] = prop1
            self.setncattr(VAR_NS_PROP1, 'name', 'ID')
            self.fh.createVariable(VAR_NS_NAMES, CHAR,
                                   (DIM_NUM_NS, DIM_LEN_STRING))
            i = 1
            for (name, nodes) in nodesets.items():
                self.fh.variables[VAR_NS_NAMES][i-1,:] = adjstr(name)
                self.fh.createDimension(DIM_NUM_NOD_NS(i), len(nodes))
                self.fh.createVariable(VAR_NOD_NS(i), INT, (DIM_NUM_NOD_NS(i),))
                self.fh.variables[VAR_NOD_NS(i)][:] = nodes + 1
                i += 1

        # ------------------------------------------------------------------- #
        # -------------------------------------------- side set meta data --- #
        # ------------------------------------------------------------------- #
        if sidesets:
            nss = len(sidesets)
            self.fh.createDimension(DIM_NUM_SS, nss)
            # side set IDs - standard map
            prop1 = arange(nss, dtype=int32) + 1
            self.fh.createVariable(VAR_SS_PROP1, INT, (DIM_NUM_SS,))
            self.fh.variables[VAR_SS_PROP1][:] = prop1
            self.setncattr(VAR_SS_PROP1, 'name', 'ID')
            self.fh.createVariable(VAR_SS_NAMES, CHAR, (DIM_NUM_SS, DIM_LEN_STRING))
            i = 1
            for (name, ss) in sidesets.items():
                ss_elems = [s[0]+1 for s in ss]
                ss_sides = [s[1]+1 for s in ss]
                self.fh.variables[VAR_SS_NAMES][i-1,:] = adjstr(name)

                self.fh.createDimension(DIM_NUM_SIDE_SS(i), len(ss_sides))
                self.fh.createVariable(VAR_SIDE_SS(i), INT, (DIM_NUM_SIDE_SS(i),))
                self.fh.variables[VAR_SIDE_SS(i)][:] = ss_sides

                self.fh.createDimension(DIM_NUM_ELE_SS(i), len(ss_elems))
                self.fh.createVariable(VAR_ELE_SS(i), INT, (DIM_NUM_ELE_SS(i),))
                self.fh.variables[VAR_ELE_SS(i)][:] = ss_elems
                i += 1

    def initialize(self, nodvarnames, elevarnames):
        '''Writes the initialization parameters to the EXODUS II file'''

        numblk = self.getdim(DIM_NUM_ELEBLK)
        # ------------------------------------------------------------------- #
        # ------------------------------------------ global variable data --- #
        # ------------------------------------------------------------------- #
        self.fh.createDimension(DIM_NUM_GLO_VAR, 1)
        self.fh.createVariable(VALS_GLO_VAR, FLOAT, (DIM_TIME_STEP, ))

        # ------------------------------------------------------------------- #
        # ----------------------------------------------------- node data --- #
        # ------------------------------------------------------------------- #
        # assert 'displx' in nodvarnames
        # make sure displacement names are in right order
        if nodvarnames:
            numnodvar = len(nodvarnames)
            nodvarnames = sorted(nodvarnames, key=sortexoname)
            self.fh.createDimension(DIM_NUM_NOD_VAR, numnodvar)
            self.fh.createVariable(VAR_NAME_NOD_VAR, CHAR,
                                   (DIM_NUM_NOD_VAR, DIM_LEN_STRING))
            for (i, nodvar) in enumerate(nodvarnames):
                key = adjstr(nodvar)
                self.fh.variables[VAR_NAME_NOD_VAR][i,:] = key
                self.fh.createVariable(VALS_NOD_VAR(i+1), FLOAT,
                                       (DIM_TIME_STEP, DIM_NUM_NOD))

        # ------------------------------------------------------------------- #
        # -------------------------------------------------- element data --- #
        # ------------------------------------------------------------------- #
        if elevarnames:
            numelevar = len(elevarnames)
            elevarnames = sorted(elevarnames, key=sortexoname)
            self.fh.createDimension(DIM_NUM_ELE_VAR, numelevar)
            self.fh.createVariable(VAR_NAME_ELE_VAR, CHAR,
                                   (DIM_NUM_ELE_VAR, DIM_LEN_STRING))
            for (i, elevar) in enumerate(elevarnames):
                key = adjstr(elevar)
                self.fh.variables[VAR_NAME_ELE_VAR][i,:] = key
                for j in range(numblk):
                    d1 = DIM_NUM_EL_IN_BLK(j+1)
                    self.fh.createVariable(VALS_ELE_VAR(i+1,j+1),
                                           FLOAT, (DIM_TIME_STEP, d1))

        self.initialized = True
        return

    def snapshot(self, step):
        if not self.initialized:
            nodvarnames, elevarnames = [], []
            for fo in step.frames[0].field_outputs.values():
                if fo.position == NODE:
                    nodvarnames.extend(fo.keys)
                else:
                    if any(in1d(fo.keys, elevarnames)):
                        continue
                    elevarnames.extend(fo.keys)
            self.initialize(nodvarnames, elevarnames)

        numele = self.getdim(DIM_NUM_ELE)
        numnod = self.getdim(DIM_NUM_NOD)
        numblk = self.getdim(DIM_NUM_ELEBLK)
        numdim = self.getdim(DIM_NUM_DIM)

        for frame in step.frames:
            if not frame.converged:
                logging.warn('CANNOT WRITE UNCONVERGED FRAME')
                return
            self.putframe(frame)

    def putframe(self, frame):
        # write time value
        count = self.count
        self.fh.variables[VAR_TIME_WHOLE][count] = frame.value
        self.fh.variables[VALS_GLO_VAR][count] = frame.increment

        nodvars = self.fh.variables.get(VAR_NAME_NOD_VAR)
        if nodvars is not None:
            nodvars = stringify2(nodvars)
            for fo in frame.field_outputs.values():
                if fo.position != NODE:
                    continue
                for (k, label) in enumerate(fo.keys):
                    i = nodvars.index(label) + 1
                    self.fh.variables[VALS_NOD_VAR(i)][count] = fo.data[:,k]

        elevars = self.fh.variables.get(VAR_NAME_ELE_VAR)
        if elevars is not None:
            elevars = stringify2(elevars)
            ebs = stringify2(self.fh.variables[VAR_EB_NAMES][:])
            for (name, fo) in frame.field_outputs.items():
                if fo.position == NODE:
                    continue
                ieb = ebs.index(name[0])
                data = fo.get_data(position=ELEMENT_CENTROID)
                for (i, label) in enumerate(fo.keys):
                    j = elevars.index(label) + 1
                    self.fh.variables[VALS_ELE_VAR(j,ieb+1)][count] = data[:,i]

        self.count += 1

class EXOFileReader(EXOFile):
    mode = 'r'
    def __init__(self, filename):
        if not isfile(filename):
            raise IOError('NO SUCH FILE: {0}'.format(repr(filename)))
        self.filename = filename
        self.fh = self.open_file(filename, mode='r')
        self.read()
        if self.numnodvar or self.numelevar:
            self.read_results()

    def read(self):
        # --- Read in the mesh

        # Basic dimensioning information
        numnod = self.getdim(DIM_NUM_NOD)
        numdim = self.getdim(DIM_NUM_DIM)
        numele = self.getdim(DIM_NUM_ELE)
        numblk = self.getdim(DIM_NUM_ELEBLK)
        numns = self.getdim(DIM_NUM_NS, 0)
        numes = self.getdim(DIM_NUM_ES, 0)
        numss = self.getdim(DIM_NUM_SS, 0)
        maxnod = max([self.getdim(DIM_NUM_NOD_PER_EL(i+1)) for i in range(numblk)])

        # Node and element maps
        # maps from external to internal numbers
        if VAR_NOD_MAP(1) in self.fh.variables:
            nm = self.fh.variables[VAR_NOD_MAP(1)]
        elif VAR_NOD_NUM_MAP in self.fh.variables:
            nm = self.fh.variables[VAR_NOD_NUM_MAP]
        elif VAR_NOD_MAP('') in self.fh.variables:
            nm = self.fh.variables[VAR_NOD_MAP('')]
        else:
            nm = range(numnod)
        nodmap = dict([(xn, n) for (n, xn) in enumerate(nm)])

        # maps from external to internal numbers
        if VAR_ELE_MAP(1) in self.fh.variables:
            em = self.fh.variables[VAR_ELE_MAP(1)]
        elif VAR_ELE_NUM_MAP in self.fh.variables:
            em = self.fh.variables[VAR_ELE_NUM_MAP]
        elif VAR_ELE_MAP('') in self.fh.variables:
            em = self.fh.variables[VAR_ELE_MAP('')]
        else:
            em = range(numele)
        elemap = dict([(xe, e) for (e, xe) in enumerate(em)])
        elemap1 = array(sorted(elemap.keys(), key=lambda k: elemap[k]))

        # Coordinates
        coord = column_stack([self.fh.variables[VAR_COOR_NAME(i)][:]
                              for i in range(numdim)])

        # ------------------------------------------------------------------- #
        # ------------------------------- element blocks and connectivity --- #
        # ------------------------------------------------------------------- #
        element_blocks = []
        elemsets = {}
        blknams = stringify2(self.fh.variables[VAR_EB_NAMES][:])
        k = 0
        for ieb in range(numblk):
            name = blknams[ieb]
            blkcon = self.fh.variables[VAR_BLKCON(ieb+1)][:]-1
            ix = arange(k, k+blkcon.shape[0])
            elefam = element_family(numdim, blkcon.shape[1])
            blk = element_block(name, len(element_blocks)+1, elemap1[ix], elefam, blkcon)
            element_blocks.append(blk)
            elemsets[name] = ix
            k += ix.shape[0]

        # ------------------------------------------------------------------- #
        # ------------------------------------------element set meta data --- #
        # ------------------------------------------------------------------- #
        elemsets = {}
        if numes:
            esnames = stringify2(self.fh.variables[VAR_ES_NAMES][:])
            for ies in range(numes):
                name = esnames[ies]
                elemsets[name] = self.fh.variables[VAR_ELE_ES(ies+1)][:] - 1

        # ------------------------------------------------------------------- #
        # -------------------------------------------- node set meta data --- #
        # ------------------------------------------------------------------- #
        nodesets = {}
        if numns:
            nsnames = stringify2(self.fh.variables[VAR_NS_NAMES][:])
            for ins in range(numns):
                name = nsnames[ins]
                nodesets[name] = self.fh.variables[VAR_NOD_NS(ins+1)][:] - 1

        # ------------------------------------------------------------------- #
        # -------------------------------------------- side set meta data --- #
        # ------------------------------------------------------------------- #
        sidesets = {}
        if numss:
            ssnames = stringify2(self.fh.variables[VAR_SS_NAMES][:])
            for iss in range(numss):
                name = ssnames[iss]
                ss_elems = self.fh.variables[VAR_ELE_SS(iss+1)][:] - 1
                ss_sides = self.fh.variables[VAR_SIDE_SS(iss+1)][:] - 1
                sidesets[name] = list(zip(ss_elems, ss_sides))

        # Save information read to self
        self.nodmap = nodmap
        self.coord = coord
        self.elemap = elemap
        self.nodesets = nodesets
        self.blknams = blknams
        self.element_blocks = element_blocks
        self.sidesets = sidesets
        self.elefam = elefam
        self.elemsets = elemsets
        self.numnod = numnod
        self.numdim = numdim
        self.numele = numele
        self.numblk = numblk
        self.numns = numns
        self.numes = numes
        self.numss = numss
        self.numelevar = self.getdim(DIM_NUM_ELE_VAR,0)
        self.numnodvar = self.getdim(DIM_NUM_NOD_VAR,0)
        self.inodmap = sorted(nodmap.keys(), key=lambda k: nodmap[k])
        self.ielemap = sorted(elemap.keys(), key=lambda k: elemap[k])

    def parse_names_and_components(self, names):
        scalars, tensors, vectors = [], {}, {}
        for (i, name) in enumerate(names):
            if name.lower().endswith(('xx', 'yy', 'zz', 'xy', 'yz', 'xz')):
                # Could be a tensor quantity
                key, c = name[:-2], name[-2:]
                tensors.setdefault(key.rstrip('_'), []).append((i, c))
            elif name.lower().endswith(('x', 'y', 'z')):
                # Could be a vector quantity
                key, c = name[:-1], name[-1]
                vectors.setdefault(key.rstrip('_'), []).append((i, c))
            else:
                # Scalar
                scalars.append((i, name))
        def key1(a):
            return {'xx':0, 'yy':1, 'zz':2, 'xy':3, 'yz':4, 'xz':5}[a[1].lower()]
        def key2(a):
            return {'x':0, 'y':1, 'z':2}[a[1].lower()]
        for (key, val) in tensors.items():
            if len(val) == 1:
                i = val[0][0]
                scalars.append((i, names[i]))
            else:
                tensors[key] = sorted(val, key=key1)
        for (key, val) in vectors.items():
            if len(val) == 1:
                i = val[0][0]
                scalars.append((i, names[i]))
            else:
                vectors[key] = sorted(val, key=key2)
        return scalars, vectors, tensors

    def read_results(self):
        """Read the results from the output database. """

        if not (self.numnodvar + self.numelevar):
            raise ValueError('NO RESULTS')

        times = self.fh.variables[VAR_TIME_WHOLE][:]
        dtimes = self.fh.variables[VALS_GLO_VAR][:]
        elevarnames = self.fh.variables.get(VAR_NAME_ELE_VAR)
        if elevarnames is not None:
            elevarnames = stringify2(elevarnames)
        else:
            elevarnames = ''
        nodvarnames = self.fh.variables.get(VAR_NAME_NOD_VAR)[:]
        if nodvarnames is not None:
            nodvarnames = stringify2(nodvarnames)
        else:
            nodvarnames = ''
        numstep = len(times)

        node_labels = sorted(self.nodmap, key=lambda k: self.nodmap[k])

        scalars1, vectors1, tensors1 = self.parse_names_and_components(nodvarnames)
        scalars2, vectors2, tensors2 = self.parse_names_and_components(elevarnames)

        self.steps = StepRepository1()
        step = self.steps.Step()
        frame = step.frames[0]

        # --- REGISTER VARIABLES

        # NODE DATA
        for (i, name) in scalars1:
            frame.ScalarField(name, NODE, node_labels)
        for (name, item) in vectors1.items():
            if name == 'displ': name = 'U'
            frame.VectorField(name, NODE, node_labels, self.numdim)
        for (name, item) in tensors1.items():
            ndir, nshr = {1:(1,0), 3:(2,1), 4:(3,1), 6:(3,3)}[len(item)]
            frame.SymmetricTensorField(name, NODE, node_labels, ndir, nshr)

        # ELEMENT DATA
        for (ieb, eb) in enumerate(self.element_blocks):
            elems = [eb.elefam] * len(eb.labels)
            for (i, name) in scalars2:
                frame.ScalarField(name, ELEMENT_CENTROID, eb.labels, eb.name)
            for (name, item) in vectors2.items():
                frame.VectorField(name, ELEMENT_CENTROID,
                                  eb.labels, self.numdim, eb.name)
            for (name, item) in tensors2.items():
                ndir, nshr = {1:(1,0), 3:(2,1), 4:(3,1), 6:(3,3)}[len(item)]
                frame.SymmetricTensorField(name, ELEMENT_CENTROID, eb.labels,
                                           ndir, nshr, eleblk=eb.name,
                                           elements=elems)

        for (count, time) in enumerate(times):
            if count > 0:
                frame = step.Frame(dtimes[count])

            # NODE DATA
            for (i, name) in scalars1:
                d = self.fh.variables[VALS_NOD_VAR(i+1)][count]
                frame.field_outputs[name].add_data(d)
            for (name, item) in vectors1.items():
                if name == 'displ': name = 'U'
                d = []
                for (i, comp) in item:
                    d.append(self.fh.variables[VALS_NOD_VAR(i+1)][count])
                d = column_stack(d)
                frame.field_outputs[name].add_data(d)
            for (name, item) in tensors1.items():
                d = []
                for (i, comp) in item:
                    d.append(self.fh.variables[VALS_NOD_VAR(i+1)][count])
                d = column_stack(d)
                frame.field_outputs[name].add_data(d)

            # ELEMENT DATA
            for (ieb, eb) in enumerate(self.element_blocks):
                for (i, name) in scalars2:
                    d = self.fh.variables[VALS_ELE_VAR(i+1,ieb+1)][count]
                    frame.field_outputs[eb.name,name].add_data(d)
                for (name, item) in vectors2.items():
                    d = []
                    for (i, comp) in item:
                        d.append(self.fh.variables[VALS_ELE_VAR(i+1,ieb+1)][count])
                    d = column_stack(d)
                    frame.field_outputs[eb.name,name].add_data(d)
                for (name, item) in tensors2.items():
                    d = []
                    for (i, comp) in item:
                        d.append(self.fh.variables[VALS_ELE_VAR(i+1,ieb+1)][count])
                    d = column_stack(d)
                    ndir, nshr = {1:(1,0), 3:(2,1), 4:(3,1), 6:(3,3)}[len(item)]
                    frame.field_outputs[eb.name,name].add_data(d)

        return self.steps

    def get_elem_coord(self, xel):
        iel = self.elemap[xel]
        k = 0
        for (ieb, eb) in enumerate(self.element_blocks):
            if iel < k+eb.elecon.shape[0]:
                inodes = eb.elecon[k+iel]
                break
            k += eb.elecon.shape[0]
        else:
            raise ValueError('COULD NOT DETERMINE NODE NUMBERS')
        return self.coord[inodes]

def put_nodal_solution(filename, nodmap, elemap, coord, element_blocks, u):

    # create new file
    exo = EXOFileWriter(filename)
    fh = exo.fh

    # initialize file with parameters
    numnod, numdim = coord.shape
    numele = len(elemap)
    exo.genesis(nodmap, elemap, coord, element_blocks)
    fh.createDimension(DIM_NUM_GLO_VAR, 1)
    fh.createVariable(VALS_GLO_VAR, FLOAT, (DIM_TIME_STEP, ))

    nodvarnames = VAR_DISP_NAMES(numdim)
    numnodvar = len(nodvarnames)
    fh.createDimension(DIM_NUM_NOD_VAR, numnodvar)
    fh.createVariable(VAR_NAME_NOD_VAR, CHAR, (DIM_NUM_NOD_VAR, DIM_LEN_STRING))
    for (k, nodvar) in enumerate(nodvarnames):
        key = adjstr(nodvar)
        fh.variables[VAR_NAME_NOD_VAR][k,:] = key
        fh.createVariable(VALS_NOD_VAR(k+1), FLOAT, (DIM_TIME_STEP, DIM_NUM_NOD))

    u0 = zeros_like(u)
    fh.variables[VAR_TIME_WHOLE][0] = 0.
    fh.variables[VALS_GLO_VAR][0] = 0.
    for (k, label) in enumerate(nodvarnames):
        fh.variables[VALS_NOD_VAR(k+1)][0] = u0[:,k]

    fh.variables[VAR_TIME_WHOLE][1] = 1.
    fh.variables[VALS_GLO_VAR][1] = 1.
    for (k, label) in enumerate(nodvarnames):
        fh.variables[VALS_NOD_VAR(k+1)][1] = u[:,k]

    exo.update()
    exo.close()
