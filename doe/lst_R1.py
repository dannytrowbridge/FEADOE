import sys
import os
import re
import math
import string
import copy
import numpy as np

#http://wiki.python.org/moin/Generators

def edge_walker(el, ose) :
    sp = ose
    fst = filter(lambda f: f[0:len(sp)] == sp, el)

    print 'GEN fst NOW =', fst
    sp = fst[0][len(sp):]
    print 'sp = ', sp
    print 'ose = ', ose

    if( sp == osp ) :
        print 'BEFORE RETURN SP=OSP'
        return
    else :
        print 'BEFORE YIELD : ', fst[0]
        yield fst[0]
    pass
    
pass





#http://stackoverflow.com/questions/1665667/python-list-filtering-and-transformation

libs = ['libIce.so.33', 'libIce.so.3.3.1', 'libIce.so.32', 'libIce.so.3.2.0']
regex = re.compile('libIce.so\.([0-9]+\.[0-9]+\.[0-9]+)')

versions = [m.group(1) for m in [regex.match(lib) for lib in libs] if m]
print versions

#http://stackoverflow.com/questions/3640359/regular-expressions-search-in-list

print filter(regex.match, libs)
print filter(regex.search, libs)


# REGEX ONLY WORKS IF THESE ARE STRINGS
kl= ['(0, 0, -1, 0, 0, 1)', '(0, 0, 0, 1, 0, 0)', '(0, 1, -1, 0, 1, 1)', '(0, 1, 0, 0, 0, 0)', '(1, 0, -1, 1, 0, 1)', '(1, 0, 0, 0, 1, 0)', '(1, 1, -1, 1, 1, 1)', '(1, 1, 0, 1, 1, 0)']

mid_edge_re = re.compile( r'[(][-01],[ ]*[-01],[ ]*[0],[ ]*[-01],[ ]*[-01],[ ]*[-01][)]', re.IGNORECASE )

print 'kl=', kl

fkl = filter(mid_edge_re.search, kl)

print 'fkl=', fkl

#kt= [(0, 0, -1, 0, 0, 1), (0, 0, 0, 1, 0, 0), (0, 1, -1, 0, 1, 1), (0, 1, 0, 0, 0, 0), (1, 0, -1, 1, 0, 1), (1, 0, 0, 0, 1, 0), (1, 1, -1, 1, 1, 1), (1, 1, 0, 1, 1, 0)]
kt= [(0, 0, -1, 0, 0, 1), (0, 0, 0, 1, 0, 0), (0, 1, -1, 0, 1, 1), (0, 1, 0, 0, 0, 0), (1, 0, -1, 1, 0, 1), (1, 0, 0, 1, 1, 0), (1, 1, -1, 1, 1, 1), (1, 1, 0, 0, 1, 0)]

for t in kt :
    print t
    if( (t[2] == 0 ) and ( t[5] == 0 ) ) :
        print '>>>>>  ', t
    pass
pass


print '%' * 80

#fst = filter(lambda f: (f[2] == 0) and (f[5] == 0) , kt)
fst = filter(lambda f: f[2] == f[5] == 0 , kt)


for t in fst :
    print t
pass




osp = (0, 0, 0)

sp = osp

while True:

    fst = filter(lambda f: f[0:len(sp)] == sp, kt)

    print 'fst NOW =', fst

    sp = fst[0][len(sp):]
    print sp

    if( sp == osp ) : break
pass

print 'kt = ', kt



# FIGURE THIS OUT LATER - EGDGE WALKER
#p = edge_walker(kt, osp)
#print 'E1 =', p.next()
#print 'E2 =', p.next()
#print 'E3 =', p.next()
#print 'E4 =', p.next()


#fst = filter(lambda f: f[0:len(sp)] == sp, kt)
#print 'fst NOW =', fst
#sp = fst[0][len(sp):]
#print sp

#yield

p = 1
n = -1

U = 0
V = 1
W = 2

bv = {}
bv[-1, -1, -1, -1, -1, -1] = 'ULO'
bv[-1, -1, -1, 1, -1, -1] = 'UHI'
#ve = (p, n, n)

vs = (n, n, n) ; lve = list(vs) ; lve[U] = p ; ve = tuple(lve)



id = vs + ve

print bv[id]


id = vs + vs
print bv[id]








