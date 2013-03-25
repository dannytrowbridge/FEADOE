from pylab import *
from numpy import *
from mpl_toolkits.mplot3d import axes3d

#nx, ny, nz = (21, 3, 5)

# create 3D points
x,y,z = mgrid[0:20:21j, -1:1:3j, -2:2:5j]

print 'x=', x
print 'y=', y
print 'z=', z

xx = x.flatten()
yy = y.flatten()
zz = z.flatten()

print 'xx=', xx
print 'yy=', yy
print 'zz=', zz

# plot 3D points
fig = figure()
ax = fig.gca(projection='3d')
ax.plot(xx, yy, zz, 'o')

show()
