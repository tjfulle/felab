from numpy import *
import distmesh as dm
random.seed(190) # Always the same results

def uniform_plate(size=.05):
    fd = lambda p: dm.drectangle(p, -1, 1, -1, 1)
    fh = dm.huniform
    p, t = dm.distmesh2d(fd, fh, size, (-1, -1, 1, 1),
                         [(-1, -1),(-1, 1),(1, -1),(1, 1)])
    return p, t

def plate_with_hole(size=.05):
    fd = lambda p: dm.ddiff(dm.drectangle(p,-1,1,-1,1),
                            dm.dcircle(p,0,0,0.5))
    fh = lambda p: 0.05+0.3*dm.dcircle(p,0,0,0.5)
    p, t = dm.distmesh2d(fd, fh, size, (-1,-1,1,1),
                         [(-1,-1),(-1,1),(1,-1),(1,1)])
    return p, t

def toaba(filename, p, t):
    fh = open(filename, 'w')
    fh.write('*Node\n')
    m = len(str(len(p)))
    for (i, xy) in enumerate(p, start=1):
        x, y = [float(_) for _ in xy]
        fh.write('{0:{m}d}, {1: .18f}, {2: .18f}\n'.format(i, x, y, m=m))
    m = max([len(str(_)) for ec in t for _ in ec])
    fh.write('*Element, type=CPE3\n')
    for (i, ec) in enumerate(t, start=1):
        a, b, c = [int(_)+1 for _ in ec]
        fh.write('{0:d}, {1:{m}d}, {2:{m}d}, {3:{m}d}\n'.format(i,a,b,c,m=m))

p, t = plate_with_hole(size=.025)
toaba('plate_with_hole.inp', p, t)
