import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm

import constants
from constants import LOC, DIR, DOF, BTAGS, GTAGS, ETAGS, FTAGS, JTAGS, FACE_GRID_INDICES

class doda(object) :
    def __init__ (self, nm) :
        self.name = nm
    pass
pass


aa = doda('DO')
bb = doda('DA')

dl = []
dl.append(bb)
dl.append(aa)


s = ', '.join([ dd.name for dd in dl ])


print 's = |'+s+'|'


a = set()
a.add(4)
a.add(9)
a.add(13)

b = set()
b.add(5)
b.add(9)
b.add(3)

sl = [ a, b ]

print a.intersection(b)


c= set()
for s in sl :
    print s
    #c = c.intersection(s)
    c = c.union(s)
pass
print c


c= sl[0]
for s in sl :
    print s
    c = c.intersection(s)

pass
print c
#print sl[0]
#print sl[1]











#args = ( 'A', 'B', 'C')

## args = []
## args.append('A')
## args.append('B')
## args.append('C')

## for i, a in enumerate(args) :
##     print ' '*3, i, ')  = ', a
## pass


## q = {}

## q['DUDE',12] = 'COW'
## print q

## q['DUDE', 12] = 'NOW'
## print q

#http://docs.scipy.org/doc/numpy/reference/generated/numpy.meshgrid.html

nx, ny, nz = (21, 3, 5)
x = np.linspace(0, 20, nx)
y = np.linspace(-1, 1, ny)
#z = np.linspace(-2, 2, nz)
print 'x=', x
print 'y=', y
#print 'z=', z
#print np.meshgrid.__doc__
xv, yv = np.meshgrid(x, y)

print 'xv=', xv
print 'yv=', yv
#print 'zv=', zv

#xv, yv = np.meshgrid(x, y, sparse=True)  # make sparse output arrays
#xv
#yv




#pl = [ np.zeros(3) ] * 4
#pl[2] = np.array[ 4.5, 5.6, 6.7 ]

pl = np.zeros( (4, 3) )

pl[2] = [4.5, 5.6, 6.7 ]



print pl
m = np.array([ 0.0, 1.16666666667, -0.5 ])

pl[0] = [ 0.0, 1.5, -0.5 ]
pl[1] = [ 0.0, 1.5, 0.0 ]
pl[2] = [ 0.0, 1.0, -0.5 ]
pl[3] = [ 0.0, 1.0, 0.0 ]


vv = m - pl[0]
vv1 = pl[1] - pl[0]


print 'VV =',vv
print 'VV1 =', vv1

vv.dot(vv1)
print vv.dot(vv1)
print vv[0]*vv1[0] + vv[1]*vv1[1] + vv[2]*vv1[2]

fig = plt.figure()
ax = plt.gca(projection='3d') # RETURNS CURRENT AXIS - CREATING ONE IF NECESSARY





maxs = [-1.0 * constants.BIG_REAL, -1.0 * constants.BIG_REAL, -1.0 * constants.BIG_REAL]
mins = [ constants.BIG_REAL, constants.BIG_REAL, constants.BIG_REAL]
for p in range(4) :
    for i in range(3) :
        maxs[i] = max(pl[p][i], maxs[i])
        mins[i] = min(pl[p][i], mins[i])
    pass
    print pl[p]
    ax.scatter(pl[p][0], pl[p][1], pl[p][2], c='red', marker='o')          
pass
ax.scatter(m[0], m[1], m[2], c='green', marker='o')          


    
print 'MAXS =', maxs
print 'MINS =', mins

ax.legend()
#ax.set_xlim3d(mins[0], maxs[0])
#ax.set_ylim3d(mins[1], maxs[1])
#ax.set_zlim3d(mins[2], maxs[2])



lo = min(mins) - 0.1 * abs(min(mins))
hi = max(maxs) + 0.1 * abs(max(maxs))
print 'LO =', lo
print 'HI =', hi

ax.set_xlim3d(lo, hi)
ax.set_ylim3d(lo, hi)
ax.set_zlim3d(lo, hi)

plt.show()
