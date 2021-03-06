import sys
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import operator

import constants
from constants import LOC, DIR, DOF, BTAGS, GTAGS, ETAGS, FTAGS, JTAGS, FACE_GRID_INDICES

# PROJECTS THE S->M VECTOR ONTO THE S->E VECTOR
def project_on(sv, mv, ev) :
    dse = ev - sv
    dsm = mv - sv
    dem = mv - ev
    ldse = np.linalg.norm(dse)
    ldsm = np.linalg.norm(dsm)
    ldem = np.linalg.norm(dem)
    
    udse = 0.0
    if( ldse != 0.0 ) : udse = dse / ldse

    udsm = 0.0
    if( ldsm != 0.0 ) :  udsm = dsm / ldsm

    udem = 0.0
    if( ldem != 0.0 ) :  udem = dem / ldem
    
    dsm_on_dse = dsm.dot(udse) * udse
    return(dsm_on_dse)
pass





# THIS ASSUMES THE THREE POINTS ARE NEARLY LINEAR
# RETURNS THE DECIMAL PRECENT OF EACH HALF OF THE DIVIDE S->M->E VECTOR SERIES
# PLUS THE UNIT VECTORS |S->E|, |M->S|, and |M->E|
def distance_factors(sv, mv, ev) :
    dse = ev - sv
    dsm = mv - sv
    dem = mv - ev
    ldse = np.linalg.norm(dse)
    ldsm = np.linalg.norm(dsm)
    ldem = np.linalg.norm(dem)

    udse = 0.0
    if( ldse != 0.0 ) : udse = dse / ldse

    udsm = 0.0
    if( ldsm != 0.0 ) :  udsm = dsm / ldsm

    udem = 0.0
    if( ldem != 0.0 ) :  udem = dem / ldem
    
    tl = ldsm + ldem
    fs = 1.0 - ldsm / tl
    fe = 1.0 - ldem / tl
    #return( fs, fe, udse, udsm, udem)
    return( fs, fe )
pass



# PROJECT S->M ONTO S->E AND THE CALC THE DECIMAL PERCENT OF THE DIVIDED S->E VECTOR
# AND THE PROJECTED VECTOR
def project_and_calc_distance_factors(sv, mv, ev) :
    dse = ev - sv
    dsm = mv - sv
    dem = mv - ev
    ldse = np.linalg.norm(dse)
    ldsm = np.linalg.norm(dsm)
    ldem = np.linalg.norm(dem)

    udse = 0.0
    if( ldse != 0.0 ) : udse = dse / ldse

    udsm = 0.0
    if( ldsm != 0.0 ) :  udsm = dsm / ldsm

    udem = 0.0
    if( ldem != 0.0 ) :  udem = dem / ldem

    dsm_on_dse = dsm.dot(udse) * udse
    fs = 1.0 - np.linalg.norm(dsm_on_dse) / ldse
    fe = 1.0 - fs
    return( fs, fe, dsm_on_dse)
pass


def project_new_point_and_calc_distance_factors(sv, mv, ev) :
    print sv, '-=>', ev
    dse = ev - sv
    dsm = mv - sv
    dem = mv - ev
    ldse = np.linalg.norm(dse)
    ldsm = np.linalg.norm(dsm)
    ldem = np.linalg.norm(dem)

    udse = 0.0
    if( ldse != 0.0 ) : udse = dse / ldse

    udsm = 0.0
    if( ldsm != 0.0 ) :  udsm = dsm / ldsm

    udem = 0.0
    if( ldem != 0.0 ) :  udem = dem / ldem

    dsm_on_dse = dsm.dot(udse) * udse
    fs = 1.0 - np.linalg.norm(dsm_on_dse) / ldse
    fe = 1.0 - fs
    npp = sv + dsm_on_dse
    return( fs, fe, npp)
pass



def interp_coef_2d_field_from_corner_points(m, pl) :
    
    f01, f10, np01 = project_new_point_and_calc_distance_factors(pl[0], m, pl[1])
    f32, f23, np32 = project_new_point_and_calc_distance_factors(pl[3], m, pl[2])
    fm01, fm32 = distance_factors(np01, m, np32)

    print 'f01 = ', f01, '   f10 =', f10, '   NP01 =', np01
    print 'f32 = ', f32, '   f23 =', f23, '   NP32 =', np32
    print 'fm01 = ', fm01, '   fm32 = ', fm32
    
    #m1 = fm01 * ( f01 * disp[0] +  f10 * disp[1] ) + fm32 * ( f23 * disp[2] + f32 * disp[3] )
    #print 'M1 =', m1


    f31, f13, np31 = project_new_point_and_calc_distance_factors(pl[3], m, pl[1])
    f02, f20, np02 = project_new_point_and_calc_distance_factors(pl[0], m, pl[2])
    fm31, fm02 = distance_factors(np31, m, np02)

    print 'f31 = ', f31, '   f13 =', f13, '   NP31 =', np31
    print 'f02 = ', f02, '   f20 =', f20, '   NP02 =', np02
    print 'fm02 = ', fm02, '   fm31 = ', fm31

    #m2 = fm02 * ( f02 * disp[0] +  f20 * disp[2] ) + fm13 * ( f13 * disp[1] + f31 * disp[3] )
    
    #print 'M2 =', m2

    #m = 1.0 / 2.0 * ( m1 + m2 )

    #SOMETHING IS NOT RIGHT HERE factor is indicating coincident grid index 2 instead it should be 1

    factors = np.zeros(4)
    factors[0] = ( fm01 * f01 + fm02 * f02 ) / 2.0
    factors[1] = ( fm01 * f10 + fm31 * f13 ) / 2.0
    factors[2] = ( fm32 * f23 + fm02 * f20 ) / 2.0
    factors[3] = ( fm32 * f32 + fm31 * f31 ) / 2.0

    return(factors)
pass





p1 = np.array([0.0, 3.0, 0.0])
p2 = np.array([0.0, 0.0, 0.0])
m = np.array([4.0, 1.0, 0.0])

f1, f2, pv = project_and_calc_distance_factors(p1, m, p2)

print 'F1=', f1
print 'F2=', f2
print 'PROJ VEC = ', pv
new_point = p1 + pv
print 'NEW POINT =', new_point


ff1, ff2, nnp = project_new_point_and_calc_distance_factors(p1, m, p2)
print 'NEW POINT =', nnp

#sys.exit(0)





class doda(object) :
    def __init__ (self, nm) :
        self.name = nm
    pass
    def __repr__(self):
        s = str(self.name)
        return(s)
    pass
pass


aa = doda(4)
bb = doda(9)
cc = doda(12)
dd = doda(2)

dl = []
dl.append(aa)
dl.append(bb)
dl.append(cc)
dl.append(dd)

print dl

#s = ', '.join([ dd.name for dd in dl ])
#print 's = |'+s+'|'



dt = reduce(lambda x, y: max([x, y], key=lambda z : z.name), dl) 

#dt = reduce(lambda x, y: filter(min([x, y], key=lambda z : z.name), [x, y]), dl) 
print 'dt=', dt

#dt = max(dl, key=operator.itemgetter(1))
dt = min( dl, key=lambda x : x.name )

print 'dt.name=', dt.name
print 'dt=', dt


for n in range(12) :
    neq = n + 2
    ln = int(neq / 2.0) + neq % 2
    print 'NEQ =', neq, 'LN = ', ln
pass


s = 'A={}, B={}, C={}'

a = s.format(4, 'R', 23.6)
print a





sys.exit()


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

#pl[2] = [4.5, 5.6, 6.7 ]



#print pl
#m = np.array([ 0.0, 1.16666666667, -0.5 ])
#m = np.array([ 0.0, 1.25, -0.5 ])
m = np.array([ 0.0, 1.5, -0.5 ])


# THE FIRST POINT SHOULD BE THE POINT NEAREST THE MASTER POINT
# THE LAST POINT SHOULD THE THE ONE CATTY-CORNER (DIAGNONAL)
pl[0] = np.array([ 0.0, 1.5, -0.5 ])
pl[1] = np.array([ 0.0, 1.5, 0.0 ])
pl[2] = np.array([ 0.0, 1.0, -0.5 ])
pl[3] = np.array([ 0.0, 1.0, 0.0 ])



print 'POINTS =', pl
print 'MASTER POINT =', m

print

disp = np.zeros( (4, 3) )

disp[0] = np.array([ 1.0, 0.0, 0.0 ])
disp[1] = np.array([ 2.0, 0.0, 0.0 ])
disp[2] = np.array([ 3.0, 0.0, 0.0 ])
disp[3] = np.array([ 4.0, 0.0, 0.0 ])


print 'DISP = ', disp

#ffs = np.array([4])
ffs = interp_coef_2d_field_from_corner_points(m, pl)

print 'FACTORS = ', ffs

mdisp = ffs[0] * disp[0] + ffs[1] * disp[1] + ffs[2] * disp[2] + ffs[3] * disp[3]
print 'MDISP =', mdisp

for i in range(4) :
    print '-'*80
    m = pl[i]
    print 'I = ', i, ' ', m
    ffs = interp_coef_2d_field_from_corner_points(m, pl)

    print 'FACTORS = ', ffs

    mdisp = ffs[0] * disp[0] + ffs[1] * disp[1] + ffs[2] * disp[2] + ffs[3] * disp[3]
    print 'MDISP =', mdisp
pass







sys.exit(0)





f01, f10, np01 = project_new_point_and_calc_distance_factors(pl[0], m, pl[1])
f32, f23, np32 = project_new_point_and_calc_distance_factors(pl[3], m, pl[2])
fm01, fm32 = distance_factors(np01, m, np32)

m1 = fm01 * ( f01 * disp[0] +  f10 * disp[1] ) + fm32 * ( f23 * disp[2] + f32 * disp[3] )
print 'M1 =', m1


f13, f31, np13 = project_new_point_and_calc_distance_factors(pl[1], m, pl[3])
f02, f20, np02 = project_new_point_and_calc_distance_factors(pl[0], m, pl[2])
fm13, fm02 = distance_factors(np13, m, np02)

m2 = fm02 * ( f02 * disp[0] +  f20 * disp[2] ) + fm13 * ( f13 * disp[1] + f31 * disp[3] )

print 'M2 =', m2

m = 1.0 / 2.0 * ( m1 + m2 )


print 'M =', m



sys.exit(0)





vv = m - pl[0]
vv1 = pl[1] - pl[0]
vv2 = pl[2] - pl[0]

print 'VV =',vv
print 'VV1 =', vv1

#print vv.dot(vv1)
#print vv[0]*vv1[0] + vv[1]*vv1[1] + vv[2]*vv1[2]

vv1_uv = vv1 / np.linalg.norm(vv1)
print 'VV1 UNIT VECTOR = ', vv1_uv

vlonvv1 = vv.dot(vv1_uv)

print 'VL ON VV1 = ', vlonvv1

vonvv1 = vv1_uv * vlonvv1


print 'V ON VV1 = ', vonvv1



print '-'*70

print 'VV2 =', vv2

vv2_uv = vv2 / np.linalg.norm(vv2)
print 'VV2 UNIT VECTOR = ', vv2_uv

vlonvv2 = vv.dot(vv2_uv)

print 'VL ON VV2 = ', vlonvv2

vonvv2 = vv2_uv * vlonvv2


print 'V ON VV2 = ', vonvv2


print 'pl[0] =', pl[0]
print 'm=',m
print 'pl[1] =', pl[1]

fs, fe, u1, u2, u3 = distance_factors(pl[0], m, pl[1])



print 'FACTOR 1 = ', fs
print 'FACTOR 2 = ', fe




m=[2,3,7,8]
c=[1,5,7,4]
d=[6,9,4,5]


print 'm=', m
print 'c=', c
print 'd=', d

ll = [ c, d ]

print 'm.intersection(c) = ', set(m).intersection(set(c))
print 'm.intersection(d) = ', set(m).intersection(set(d))


print ll

nll = filter( lambda g : len(set(m).intersection(set(g))) != 0, ll)

print nll

sys.exit(0)














# PLOT STUFF


fig = plt.figure()
ax = plt.gca(projection='3d') # RETURNS CURRENT AXIS - CREATING ONE IF NECESSARY


maxs = np.array([-1.0 * constants.BIG_REAL, -1.0 * constants.BIG_REAL, -1.0 * constants.BIG_REAL])
mins = np.array([ constants.BIG_REAL, constants.BIG_REAL, constants.BIG_REAL] )
for p in range(4) :
    pp = p + 1
    if( pp > 3 ) : pp = 0
    for i in range(3) :
        maxs[i] = max(pl[p][i], maxs[i])
        mins[i] = min(pl[p][i], mins[i])
    pass
    #print pl[p]

    x = [ pl[p][0], pl[pp][0] ]
    y = [ pl[p][1], pl[pp][1] ]
    z = [ pl[p][2], pl[pp][2] ]
    ax.plot(x, y, z, zdir='z', color='black')
 
    ax.scatter(pl[p][0], pl[p][1], pl[p][2], c='red', marker='o')          
pass
ax.scatter(m[0], m[1], m[2], c='green', marker='o')          


    
#print 'MAXS =', maxs
#print 'MINS =', mins

ax.legend()
#ax.set_xlim3d(mins[0], maxs[0])
#ax.set_ylim3d(mins[1], maxs[1])
#ax.set_zlim3d(mins[2], maxs[2])



lo = min(mins) - 0.1 * abs(min(mins))
hi = max(maxs) + 0.1 * abs(max(maxs))

vlo = mins - 0.1 * abs(mins)
vhi = maxs + 0.1 * abs(maxs)
#print 'LO =', lo
#print 'HI =', hi
#print 'VLO =', vlo
#print 'VHI =', vhi

#ax.set_xlim3d(lo, hi)
#ax.set_ylim3d(lo, hi)
#ax.set_zlim3d(lo, hi)
ax.set_xlim3d(vlo[0], vhi[0])
ax.set_ylim3d(vlo[1], vhi[1])
ax.set_zlim3d(vlo[2], vhi[2])

plt.show()
