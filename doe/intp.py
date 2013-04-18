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
    #print 'VALS =', vals
    #print 'pt #', i, ' = ', pt[i]
    #gvv = griddata(pt, vals, (3.5, 4.5), method='nearest')
    #gvv = griddata(pt, vals, (3.5, 4.5), method='linear')
    gvv = griddata(pt, vals, (3.5, 4.5), method='cubic')
    #print 'GVV=', gvv
    sum = sum + gvv
    gv.append(gvv)
pass

#print 'GV=', gv
#print 'sum =', sum 


v = [ 1, 2, 3 ]

print type(v)
print isinstance(v, list)

n=-1
p=1

e = {}

e[n,n,n,p,n,n] = '1_to_2'
e[p,n,n,p,p,n] = '2_to_3'
e[n,n,n,n,p,n] = '1_to_4'
e[n,p,n,p,p,n] = '4_to_1'

pto = [ (None, None, None),
        (n, n, n),
        (p, n, n),
        (p, p, n),
        (n, p, n),
        (n, n, p),
        (p, n, p),
        (p, p, p),
        (n, p, p) ]


tt23 = ( pto[2] + pto[3] )
print tt23
print 'E23 =', e[tt23]
print 'E23 =', e[pto[2] + pto[3]]
print 'E23 =', e[ ( pto[2] + pto[3] ) ]


ll =  ()
print '_' * 80
print 'll =', ll, ' type =', type(ll)

if( not isinstance(ll, list) ) :
    print ' Change ll to a list:', ll
    if( isinstance(ll, tuple) ) :
        ll = list(ll)
    else :
        ll = [ ll ]
    pass
    print '  NOW ll :', ll
pass

print 'll =', ll
