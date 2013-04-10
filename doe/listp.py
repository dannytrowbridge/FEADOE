
import constants
from constants import LOC, DIR, DOF, BTAGS, GTAGS, ETAGS, FTAGS, JTAGS, FACE_GRID_INDICES



            
p = 1
n = -1


til = [ [ 1, (n,n,n), [ [ FTAGS.MIN_U, [ 3, 2 ] ],
                        [ FTAGS.MIN_V, [ 0, 3 ] ],
                        [ FTAGS.MIN_W, [ 0, 3 ] ],
                        ] ],
        [ 2, (p,n,n), [ [ FTAGS.MAX_U, [ 0, 3 ] ],
                        [ FTAGS.MIN_V, [ 1, 0 ] ],
                        [ FTAGS.MIN_W, [ 1, 0 ] ],
                        ] ],
        [ 3, (p,p,n), [ [ FTAGS.MAX_U, [ 3, 2 ] ],
                        [ FTAGS.MAX_V, [ 0, 3 ] ],
                        [ FTAGS.MIN_W, [ 2, 1 ] ],
                        ] ],
        [ 4, (n,p,n), [ [ FTAGS.MIN_U, [ 0, 3 ] ],
                        [ FTAGS.MAX_V, [ 3, 2 ] ],
                        [ FTAGS.MIN_W, [ 3, 2 ] ],
                        ] ],
        [ 5, (n,n,p), [ [ FTAGS.MIN_U, [ 2, 1 ] ],
                        [ FTAGS.MIN_V, [ 1, 0 ] ],
                        [ FTAGS.MAX_W, [ 0, 3 ] ],
                        ] ],
        [ 6, (p,n,p), [ [ FTAGS.MAX_U, [ 1, 0 ] ],
                        [ FTAGS.MIN_V, [ 2, 1 ] ],
                        [ FTAGS.MAX_W, [ 3, 2 ] ],
                        ] ],
        [ 7, (p,p,p), [ [ FTAGS.MAX_U, [ 2, 1 ] ],
                        [ FTAGS.MAX_V, [ 1, 0 ] ],
                        [ FTAGS.MAX_W, [ 2, 1 ] ],
                        ] ],
        [ 8, (n,p,p), [ [ FTAGS.MIN_U, [ 1, 0 ] ],
                        [ FTAGS.MAX_V, [ 2, 1 ] ],
                        [ FTAGS.MAX_W, [ 1, 0 ] ],
                        ] ]
        ]


## for n, pi, fiv in til :
##     print 'n=', n
##     print 'pi=', pi
##     print 'fiv=', fiv
##     for ft, il in fiv:
##         print 'FACE TAG =', ft,
##         print '  START POINT EDGE INDEX =', il[0], 
##         print '  END POINT EDGE INDEX =', il[1]
##     pass
        
## pass


kl = [ (-1, 1, -1, 1, 1, -1),
       (-1, 1, 1, 1, 1, 1),
       (-1, -1, -1, 1, -1, -1),
       (-1, 1, -1, -1, 1, 1),
       (-1, -1, 1, -1, 1, 1),
       (1, -1, 1, 1, 1, 1),
       (1, -1, -1, 1, -1, 1),
       (-1, -1, -1, -1, -1, 1),
       (1, 1, -1, 1, 1, 1),
       (-1, -1, -1, -1, 1, -1),
       (-1, -1, 1, 1, -1, 1),
       (1, -1, -1, 1, 1, -1)
       ]


tv = ( 1, -1, 1)

for k in kl :
    print '|', k[0:3], '|', tv, '|',  '   FIRST EQ? =', (k[0:3]==tv)
    print '|', k[3:], '|', tv, '|',  '   SECOND EQ? =', (k[3:]==tv)
pass


fkl1 = filter( ( lambda k: k[0:3] == tv ), kl)

print fkl1

fkl1 = filter( ( lambda k: k[3:] == tv ), kl)

print fkl1
