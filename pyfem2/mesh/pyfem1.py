import re
from numpy import array

def ReadInput(filename, output=None, ptabs=None):

    string = open(filename, 'r').read()

    # FIND NODES
    match = re.search(r'(?ims)<nodes>(?P<n>.*)?<\/nodes>', string)
    if not match:
        raise ValueError('NO NODES FOUND')
    labels, coords = [], []
    for line in match.group('n').split('\n'):
        for item in line.split(';'):
            item = item.split()
            if not item:
                continue
            labels.append(int(item[0]))
            coords.append([float(x) for x in item[1:]])
    labels = array(labels)
    offset = 1 if labels.min() == 0 else 0
    nmap = dict([(label, label+offset) for label in labels])
    nodtab = [[label+offset]+coords[i] for (i, label) in enumerate(labels)]

    # FIND ELEMENTS
    match = re.search(r'(?ims)<elements>(?P<e>.*)?<\/elements>', string)
    if not match:
        raise ValueError('NO ELEMENTS FOUND')
    eletab = []
    for line in match.group('e').split('\n'):
        for item in line.split(';'):
            item = item.split()
            if not item:
                continue
            label = int(item[0])
            eletyp = item[1]
            connect = array([int(x) for x in item[2:]])
            if len(connect) == 8:
                # REORDER FOR PYFEM2
                connect = connect[[0,2,4,6,1,3,5,7]]
            elif len(connect) == 6:
                # REORDER FOR PYFEM2
                connect = connect[[0,2,4,1,3,5]]
            eletab.append([label] + [nmap[x] for x in connect])

    # BOUNDARY CONDITIONS
    rx = r'(?ims)<nodeconstraints>(?P<nc>.*)?<\/nodeconstraints>'
    match = re.search(rx, string)
    if match:
        bcs = {}
        for line in match.group('nc').split('\n'):
            for item in line.split(';'):
                if not item.split():
                    continue
                item = item.split('=')
                a, n = re.split(r'[\[\]]', item[0].strip())[:2]
                bcs[(nmap[int(n)], a.upper())] = float(item[1])

    rx = r'(?ims)<externalforces>(?P<cf>.*)?<\/externalforces>'
    match = re.search(rx, string)
    if match:
        cfs = {}
        for line in match.group('cf').split('\n'):
            for item in line.split(';'):
                if not item.split():
                    continue
                item = item.split('=')
                a, n = re.split(r'[\[\]]', item[0].strip())[:2]
                cfs[(nmap[int(n)], a.upper())] = float(item[1])

    if output is None:
        output = splitext(filename) + '.inp'

    fh = open(output, 'w')
    fh.write('*Node, nset=NAll\n')
    for node in nodtab:
        fh.write(' {0:d}, '.format(node[0]))
        fh.write(', '.join('{0: .8f}'.format(float(x)) for x in node[1:]))
        fh.write('\n')

    fh.write('*Element, type=None, elset=EAll\n')
    for elem in eletab:
        fh.write(' ' + ', '.join('{0:d}'.format(int(x)) for x in elem))
        fh.write('\n')

    fh.write('*Solid Section, elset=EAll, material=None\n')

    fh.write('*Step, name=Step-1\n*Static\n')

    if bcs:
        fh.write('*Boundary\n')
        for key in sorted(bcs.keys(), key=lambda x: x[0]):
            line = ' {0:d}'.format(key[0])
            dof = {'U': ', 1, 1', 'V': ', 2, 2'}[key[1]]
            line += dof + ', {0:f}'.format(bcs[key])
            fh.write(line + '\n')

    if cfs:
        fh.write('*CLoad\n')
        for key in sorted(cfs.keys(), key=lambda x: x[0]):
            line = ' {0:d}'.format(key[0])
            dof = {'U': ', X', 'V': ', Y'}[key[1]]
            line += dof + ', {0:f}'.format(cfs[key])
            fh.write(line + '\n')

    fh.write('*End Step\n')

    fh.close()

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('source')
    p.add_argument('-o')
    p.add_argument('-p', action='store_true')
    args = p.parse_args()
    ReadInput(args.source, args.o, args.p)
