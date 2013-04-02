import sys
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import operator
import scipy 
from scipy.interpolate import griddata

pt = np.zeros((8,2))
pt[0] = np.array( [2.0, 3.0] )
pt[1] = np.array( [6.0, 3.0] )
pt[2] = np.array( [6.0, 8.0] )
pt[3] = np.array( [2.0, 8.0] )
pt[4] = np.array( [4.0, 3.0] )
pt[5] = np.array( [6.0, 5.5] )
pt[6] = np.array( [4.0, 8.0] )
pt[7] = np.array( [2.0, 5.5] )


vals = np.zeros(8)
vals[7] = 1.0

print pt

# griddata didn't like it when I had 3d points and all the z values were zero
gv = []

sum = 0.0
for i in range(8) :
    vals = np.zeros(8)
    vals[i] = 1.0
    print vals
    print pt[i]
    #gvv = griddata(pt, vals, (3.5, 4.5), method='nearest')
    #gvv = griddata(pt, vals, (3.5, 4.5), method='linear')
    gvv = griddata(pt, vals, (3.5, 4.5), method='cubic')
    print gvv
    sum = sum + gvv
    gv.append(gvv)
pass

print gv
print 'sum =', sum 
