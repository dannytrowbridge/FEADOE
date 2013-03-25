



import sys
import os
import re
import math
import string
import copy
import numpy as np
from collections import OrderedDict

##########################################################################################################

 
def area_of_triangle(p1, p2, p3):
    return np.linalg.norm(np.cross((p2.v - p1.v), (p3.v - p1.v)))/2.0
pass

##########################################################################################################

def unit_normal(p0, pu, pv):
    # NORMAL FROM CROSS PRODUCT
    pw = np.cross( (pu.v - p0.v), (pv.v - p0.v) )
    # UNIT VECTOR NORMAL
    npw = pw / np.linalg.norm(pw)
    return(npw)
pass


##########################################################################################################

class point(object) :

    # THESE ARE PARAMS ADDED TO EVERY POINT - WITH SOME DEFAULT VALUES
    # CLASS-LEVEL VARIABLE
    standard_params = [
        [ 'TEMP', 0.0 ],
        [ 'PRES', 0.0 ],
        [ 'THICK', 1.0 ], 
        [ 'VOL_FRAC', 1.0 ]
        ]

    #----------------------------------------------------------------------------------------------------

    def __init__ (self, x = 0.0, y = 0.0, z = 0.0, pp = {}) :
        self.v = np.array([x, y, z])
        self.params = OrderedDict(point.standard_params)
        #self.params = OrderedDict(self.__class__.standard_params)
        for k, v  in pp.items() :
            self.params[k] = v
        pass
    pass

    #----------------------------------------------------------------------------------------------------

    def __del__ (self) :
        del self.v
        del self.params
    pass

    #----------------------------------------------------------------------------------------------------

    def __repr__(self):
        b2b_brace_re = re.compile( r'[}][ ]*[{]', re.IGNORECASE )
        s = 'POINT: X={} Y={} Z={} ' + '{}={} ' * len(self.params)
        s =  b2b_brace_re.sub(r'}, {', s)
        vl = []
        for k, v in self.params.items() :
            vl.append(k)
            vl.append(v)
        pass
        s = s.format(self.v[0], self.v[1], self.v[2], *vl)
        return(s)
    pass
    
    #----------------------------------------------------------------------------------------------------

    def set_param(self, pname, val ) :
        self.params[pname] = val
    pass

    #----------------------------------------------------------------------------------------------------

    def get_param(self, pname, default=None) :
        val = self.params.get(pname, default)
        return(val)
    pass

    #----------------------------------------------------------------------------------------------------

    def has_param(self, pname) :
        rc = pname in self.params
        return(rc)
    pass

pass



##########################################################################################################
# http://stackoverflow.com/questions/682504/what-is-a-clean-pythonic-way-to-have-multiple-constructors-in-python

class vector(object) :

#class Cheese():
#    def __init__(self, *args, **kwargs):
#        #args -- tuple of anonymous arguments
#        #kwargs -- dictionary of named arguments
#        self.num_holes = kwargs.get('num_holes',random_holes())
        
#    def __init__ (self, x1, y1, z1, x2, y2, z2) :
    def __init__ (self, p1, p2 = None) :

        # ONLY ONE POINT SPECIFIED THEN ASSUME THIS IS A LOCATION VECTOR
        #                    - START POINT ASSUMED TO BE ORIGIN
        if( p2 is None ) : 
            self.ep = copy.deepcopy(p1)
            self.sp = point(0.0, 0.0, 0.0)

            # INITIALIZE PARAMETERS
            for k, v in self.ep.params.items() :
                # CHECK FOR STRING PLUS THE 4 PYTHON NUMERIC TYPES - OTHERWISE INITIALIZE TO None
                if( isinstance(v, str) ) :
                    self.sp.set_param(k, '')
                elif( isinstance(v, int) ) :
                    self.sp.set_param(k, 0)
                elif( isinstance(v, long) ) :
                    self.sp.set_param(k, 0L)
                elif( isinstance(v, float) ) :
                    self.sp.set_param(k, 0.0)
                elif( isinstance(v, complex) ) :
                    self.sp.set_param(k, 0.0 + 0.0J)
                else :
                    #self.sp.set_param(k], v)
                    self.sp.set_param(k, None)
                pass
            pass
        else :  
            self.sp = copy.deepcopy(p1)
            self.ep = copy.deepcopy(p2)
        pass
        #print 'START POINT:', self.sp.v
        #print 'END POINT:', self.ep.v
        self.dv = np.subtract(self.ep.v, self.sp.v)
        self.len = np.linalg.norm(self.dv)
        
        self.uv = np.array( [0.0, 0.0, 0.0] )
        if( self.len > 0.0 ) :
            self.uv = np.array( [ self.dv[0] / self.len,
                                  self.dv[1] / self.len,
                                  self.dv[2] / self.len ] )
        pass
        #print 'dv =', self.dv
        #print 'LEN = ', self.len
        #print 'uv =', self.uv, np.linalg.norm(self.uv)
        
    pass

    #----------------------------------------------------------------------------------------------------

    def __repr__(self):
        s = 'VECTOR: \n FROM: ' + str(self.sp) + '\n TO:   ' + str(self.ep) 
        s += '\n  DV = ' + str(self.dv) + '   LEN = ' + str(self.len) + '   UNIT_VEC = ' + str(self.uv)
        return(s)
    pass
    
    #----------------------------------------------------------------------------------------------------
    
    def vec_cross(self, v) : 
        npw = np.cross( self.dv, v.dv ) 
        return(npw)
    pass

    #----------------------------------------------------------------------------------------------------
    
    def vec_normal(self) :
        npw = self.dv / np.linalg.norm(self.dv)
        return(npw)
    pass

    #----------------------------------------------------------------------------------------------------

    def vec_length(self) :
        nps = np.linalg.norm(self.dv)
        return(nps)
    pass

    #----------------------------------------------------------------------------------------------------

    def scale(self, sf) :
        v = self.sp.v + self.dv * sf
        rd = OrderedDict()
        rd['X'] = v[0]
        rd['Y'] = v[1]
        rd['Z'] = v[2]
        #print 'SCALE VECTOR ', 
        #print 'ORIGINAL DV = ', self.dv
        #print 'SCALED DV = ', self.dv * sf
        #print 'NEW V = ', v
        for sk, sv in self.sp.params.items() :
            if( sv is None ) : continue
                
            if( self.sp.has_param(sk) ) :
                ev = self.ep.get_param(sk)
                if( ev is None ) : continue
                
                # http://stackoverflow.com/questions/432842/how-do-you-get-the-logical-xor-of-two-variables-in-python
                # XOR: bool(a) != bool(b)
                # IF ONE OR THE OTHER IS A STRING - BUT NOT BOTH - CONVERT THE NON-STRING TO A STRING
                #   AND HANDLE AS DISCRETE VALUES
                if( isinstance(sv, str) != isinstance(ev, str) ) :
                    if( not isinstance(sv, str) ) :
                        sv = str(sv)
                    pass
                    if( not isinstance(ev, str) ) :
                        ev = str(ev)
                    pass
                pass

                # STRING PARAMETERS...
                if( isinstance(sv, str) and isinstance(ev, str) ) :
                    rd[sk] = sv
                    if( sf > 0.5 ) : # DISCRETE CHANGE OF STRING VALUE AT PARAMETERIC MID POINT
                        rd[sk] = ev
                    pass
                    continue
                pass
                                
                dv = ev - sv
                vv = sv + dv * sf
                rd[sk] = vv
            pass
        pass
        return(rd)
    pass

    #----------------------------------------------------------------------------------------------------
    # SPECIFY A SCALE FACTOR - RETURNS A NEW POINT
    def point_from_scale(self, sf) :
        rs = self.scale(sf)
        rest = OrderedDict( [(k, rs[k]) for k in rs.keys()[3:] ] )
        p = point(rs['X'], rs['Y'], rs['Z'], rest)
        return(p)
    pass

    #----------------------------------------------------------------------------------------------------

    # SPECIFY A LENGTH - RETURNS A NEW POINT - DEFAULTS TO A UNIT VECTOR
    def point_from_length(self, new_len = 1.0) :
        sf = 1.0
        if( self.len != 0.0 ) :
            sf = new_len / self.len
        pass
        p =  self.point_from_scale(sf)
        return(p)
    pass

pass

##########################################################################################################
import sys
import os
import re
import math
import string
import copy
import numpy as np
from collections import OrderedDict

##########################################################################################################

 
def area_of_triangle(p1, p2, p3):
    return np.linalg.norm(np.cross((p2.v - p1.v), (p3.v - p1.v)))/2.0
pass

##########################################################################################################

def unit_normal(p0, pu, pv):
    # NORMAL FROM CROSS PRODUCT
    pw = np.cross( (pu.v - p0.v), (pv.v - p0.v) )
    # UNIT VECTOR NORMAL
    npw = pw / np.linalg.norm(pw)
    return(npw)
pass


##########################################################################################################

class point(object) :

    # THESE ARE PARAMS ADDED TO EVERY POINT - WITH SOME DEFAULT VALUES
    # CLASS-LEVEL VARIABLE
    standard_params = [
        [ 'TEMP', 0.0 ],
        [ 'PRES', 0.0 ],
        [ 'THICK', 1.0 ], 
        [ 'VOL_FRAC', 1.0 ]
        ]

    #----------------------------------------------------------------------------------------------------

    def __init__ (self, x = 0.0, y = 0.0, z = 0.0, pp = {}) :
        self.v = np.array([x, y, z])
        self.params = OrderedDict(point.standard_params)
        #self.params = OrderedDict(self.__class__.standard_params)
        for k, v  in pp.items() :
            self.params[k] = v
        pass
    pass

    #----------------------------------------------------------------------------------------------------

    def __del__ (self) :
        del self.v
        del self.params
    pass

    #----------------------------------------------------------------------------------------------------

    def __repr__(self):
        b2b_brace_re = re.compile( r'[}][ ]*[{]', re.IGNORECASE )
        s = 'POINT: X={} Y={} Z={} ' + '{}={} ' * len(self.params)
        s =  b2b_brace_re.sub(r'}, {', s)
        vl = []
        for k, v in self.params.items() :
            vl.append(k)
            vl.append(v)
        pass
        s = s.format(self.v[0], self.v[1], self.v[2], *vl)
        return(s)
    pass
    
    #----------------------------------------------------------------------------------------------------

    def set_param(self, pname, val ) :
        self.params[pname] = val
    pass

    #----------------------------------------------------------------------------------------------------

    def get_param(self, pname, default=None) :
        val = self.params.get(pname, default)
        return(val)
    pass

    #----------------------------------------------------------------------------------------------------

    def has_param(self, pname) :
        rc = pname in self.params
        return(rc)
    pass

pass



##########################################################################################################
# http://stackoverflow.com/questions/682504/what-is-a-clean-pythonic-way-to-have-multiple-constructors-in-python

class vector(object) :

#class Cheese():
#    def __init__(self, *args, **kwargs):
#        #args -- tuple of anonymous arguments
#        #kwargs -- dictionary of named arguments
#        self.num_holes = kwargs.get('num_holes',random_holes())
        
#    def __init__ (self, x1, y1, z1, x2, y2, z2) :
    def __init__ (self, p1, p2 = None) :

        # ONLY ONE POINT SPECIFIED THEN ASSUME THIS IS A LOCATION VECTOR
        #                    - START POINT ASSUMED TO BE ORIGIN
        if( p2 is None ) : 
            self.ep = copy.deepcopy(p1)
            self.sp = point(0.0, 0.0, 0.0)

            # INITIALIZE PARAMETERS
            for k, v in self.ep.params.items() :
                # CHECK FOR STRING PLUS THE 4 PYTHON NUMERIC TYPES - OTHERWISE INITIALIZE TO None
                if( isinstance(v, str) ) :
                    self.sp.set_param(k, '')
                elif( isinstance(v, int) ) :
                    self.sp.set_param(k, 0)
                elif( isinstance(v, long) ) :
                    self.sp.set_param(k, 0L)
                elif( isinstance(v, float) ) :
                    self.sp.set_param(k, 0.0)
                elif( isinstance(v, complex) ) :
                    self.sp.set_param(k, 0.0 + 0.0J)
                else :
                    #self.sp.set_param(k], v)
                    self.sp.set_param(k, None)
                pass
            pass
        else :  
            self.sp = copy.deepcopy(p1)
            self.ep = copy.deepcopy(p2)
        pass
        #print 'START POINT:', self.sp.v
        #print 'END POINT:', self.ep.v
        self.dv = np.subtract(self.ep.v, self.sp.v)
        self.len = np.linalg.norm(self.dv)
        
        self.uv = np.array( [0.0, 0.0, 0.0] )
        if( self.len > 0.0 ) :
            self.uv = np.array( [ self.dv[0] / self.len,
                                  self.dv[1] / self.len,
                                  self.dv[2] / self.len ] )
        pass
        #print 'dv =', self.dv
        #print 'LEN = ', self.len
        #print 'uv =', self.uv, np.linalg.norm(self.uv)
        
    pass

    #----------------------------------------------------------------------------------------------------

    def __repr__(self):
        s = 'VECTOR: \n FROM: ' + str(self.sp) + '\n TO:   ' + str(self.ep) 
        s += '\n  DV = ' + str(self.dv) + '   LEN = ' + str(self.len) + '   UNIT_VEC = ' + str(self.uv)
        return(s)
    pass
    
    #----------------------------------------------------------------------------------------------------
    
    def vec_cross(self, v) : 
        npw = np.cross( self.dv, v.dv ) 
        return(npw)
    pass

    #----------------------------------------------------------------------------------------------------
    
    def vec_normal(self) :
        npw = self.dv / np.linalg.norm(self.dv)
        return(npw)
    pass

    #----------------------------------------------------------------------------------------------------

    def vec_length(self) :
        nps = np.linalg.norm(self.dv)
        return(nps)
    pass

    #----------------------------------------------------------------------------------------------------

    def scale(self, sf) :
        v = self.sp.v + self.dv * sf
        rd = OrderedDict()
        rd['X'] = v[0]
        rd['Y'] = v[1]
        rd['Z'] = v[2]
        #print 'SCALE VECTOR ', 
        #print 'ORIGINAL DV = ', self.dv
        #print 'SCALED DV = ', self.dv * sf
        #print 'NEW V = ', v
        for sk, sv in self.sp.params.items() :
            if( sv is None ) : continue
                
            if( self.sp.has_param(sk) ) :
                ev = self.ep.get_param(sk)
                if( ev is None ) : continue
                
                # http://stackoverflow.com/questions/432842/how-do-you-get-the-logical-xor-of-two-variables-in-python
                # XOR: bool(a) != bool(b)
                # IF ONE OR THE OTHER IS A STRING - BUT NOT BOTH - CONVERT THE NON-STRING TO A STRING
                #   AND HANDLE AS DISCRETE VALUES
                if( isinstance(sv, str) != isinstance(ev, str) ) :
                    if( not isinstance(sv, str) ) :
                        sv = str(sv)
                    pass
                    if( not isinstance(ev, str) ) :
                        ev = str(ev)
                    pass
                pass

                # STRING PARAMETERS...
                if( isinstance(sv, str) and isinstance(ev, str) ) :
                    rd[sk] = sv
                    if( sf > 0.5 ) : # DISCRETE CHANGE OF STRING VALUE AT PARAMETERIC MID POINT
                        rd[sk] = ev
                    pass
                    continue
                pass
                                
                dv = ev - sv
                vv = sv + dv * sf
                rd[sk] = vv
            pass
        pass
        return(rd)
    pass

    #----------------------------------------------------------------------------------------------------
    # SPECIFY A SCALE FACTOR - RETURNS A NEW POINT
    def point_from_scale(self, sf) :
        rs = self.scale(sf)
        rest = OrderedDict( [(k, rs[k]) for k in rs.keys()[3:] ] )
        p = point(rs['X'], rs['Y'], rs['Z'], rest)
        return(p)
    pass

    #----------------------------------------------------------------------------------------------------

    # SPECIFY A LENGTH - RETURNS A NEW POINT - DEFAULTS TO A UNIT VECTOR
    def point_from_length(self, new_len = 1.0) :
        sf = 1.0
        if( self.len != 0.0 ) :
            sf = new_len / self.len
        pass
        p =  self.point_from_scale(sf)
        return(p)
    pass

pass

##########################################################################################################

def check(*el) :
    edges = [] 
    for i, edg in enumerate(el) :
        edges.append(edg)
    pass

    print edges
    
pass


##########################################################################################################



if __name__ == '__main__' :
    p1 = point(0.0, 0.0, 0.0, {'T': 20, 'VF': 1.0 } )
    p2 = point(10.0, 0.0, 0.0, {'T': 40, 'VF': 0.5 } )
    p3 = point(10.0, 10.0, 0.0, {'T': 20, 'VF': 1.0 } )
    p4 = point(0.0, 10.0, 0.0, {'T': 20, 'VF': 1.0 } )


    pl = ( p1, p2, p3, p4 )

    check(pl)
    
pass
