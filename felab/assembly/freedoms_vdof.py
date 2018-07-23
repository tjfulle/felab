from numpy import zeros, array

from ..constants import *
from ..utilities import *


def node_freedom_association_table(numnod, elements, disp=0):
    # NODE FREEDOM ASSOCIATION TABLE
    active_dof = [None] * MDOF
    nodfat = zeros((numnod, MDOF), dtype=int)
    for el in elements:
        for (i, node) in enumerate(el.inodes):
            nfs = el.signature[i]
            nf = [max(nfs[j], nodfat[node,j]) for j in range(MDOF)]
            nodfat[node] = nf
            for (j, k) in enumerate(nfs):
                if k:
                    active_dof[j] = j
    active_dof = array([x for x in active_dof if x is not None])
    if not disp:
        return nodfat
    return nodfat, active_dof


def total_degrees_of_freedom(nodfat):
    # TOTAL NUMBER OF DEGREES OF FREEDOM
    return sum(count_digits(p) for p in nodfat)


def node_freedom_map_table(nodfat, disp=0):
    # NODE FREEDOM MAP TABLE
    nodfmt = zeros(nodfat.shape[0], dtype=int)
    dofmap = {}
    dof = 0
    nodfmt = [0]
    for i in range(nodfat.shape[0]):
        for (j, k) in enumerate(nodfat[i]):
            if not k: continue
            dofmap[i,j] = dof
            dof += 1
        nodfmt.append(dof)
    nodfmt = array(nodfmt[:-1], dtype=int)
    if not disp:
        return nodfmt
    return nodfmt, dofmap


def element_freedom_table(nodfat, nodfmt, elements):
    eftab = []
    for el in elements:
        eft = zeros(sum(count_digits(nfs) for nfs in el.signature), dtype=int)
        k, count = 0, 0
        for (i, inode) in enumerate(el.inodes):
            ix, sx = 0, zeros(MDOF, dtype=int)
            nfs1 = nodfat[inode]
            for j in range(MDOF):
                if nfs1[j] > 0:
                    sx[j] = ix
                    ix += 1
            nfs2 = el.signature[i]
            for j in range(MDOF):
                if nfs2[j] > 0:
                    if nfs1[j] > 0:
                        eft[k] = nodfmt[inode] + sx[j]
                        count += 1
                    k += 1
        if all(eft==0.):
            raise UserInputError('ZERO ENTRY IN EFTAB FOR '
                                 'ELEMENT {0}'.format(el.label))
        eftab.append(eft)
    return eftab
