#http://stackoverflow.com/questions/5124126/python-scipy-interpolation-map-coordinates#

import numpy
from scipy import interpolate
x = numpy.array([0.0, 0.60, 1.0])
y = numpy.array([0.0, 0.25, 0.80, 1.0])
z = numpy.array([ 
   [ 1.4 ,  6.5 ,  1.5 ,  1.8 ],
   [ 8.9 ,  7.3 ,  1.1 ,  1.09],
   [ 4.5 ,  9.2 ,  1.8 ,  1.2 ]])
# you have to set kx and ky small for this small example dataset
# 3 is more usual and is the default
# s=0 will ensure this interpolates.  s>0 will smooth the data
# you can also specify a bounding box outside the data limits
# if you want to extrapolate
sp = interpolate.RectBivariateSpline(x, y, z, kx=2, ky=2, s=0)

sp([0.60], [0.25])  # array([[ 7.3]])
sp([0.25], [0.60])  # array([[ 2.66427408]])
