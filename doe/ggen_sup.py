import sys
import os
import re
import math
import string
import copy
import numpy as np
from collections import OrderedDict
import constants
from constants import LOC, DIR, DOF, BTAGS, GTAGS, ETAGS, FTAGS, JTAGS, FACE_GRID_INDICES

import matplotlib.pyplot as plt

#from gfea_model import mesh

##########################################################################################################
# http://wiki.python.org/moin/Generators
# http://stackoverflow.com/questions/231767/the-python-yield-keyword-explained
# http://stackoverflow.com/questions/1665667/python-list-filtering-and-transformation
#
# EDGE WALKER GENERATOR - LOOP THROUGH LIST OF INDEXES UNTIL IT COMES BACK TO WHERE IT STARTED
#   Looks for OSE in the first half of an index
#     and then loops until OSE if found in the last half of an index
#     Good for walking around the edges of a face.
##########################################################################################################
# EXAMPLE USAGE:
# kt= [(0, 0, -1, 0, 0, 1), (0, 0, 0, 1, 0, 0), (0, 1, -1, 0, 1, 1), (0, 1, 0, 0, 0, 0), (1, 0, -1, 1, 0, 1), (1, 0, 0, 1, 1, 0), (1, 1, -1, 1, 1, 1), (1, 1, 0, 0, 1, 0)]
# osp = (0, 0, 0)
## p = edge_walker(kt, osp)
## cnt = 0
## for pp in p :
##     cnt += 1
##     print 'E' + str(cnt) + ' =', pp
## pass
#
# YEILDS...
# E1 = (0, 0, 0, 1, 0, 0)
# E2 = (1, 0, 0, 1, 1, 0)
# E3 = (1, 1, 0, 0, 1, 0)
# E4 = (0, 1, 0, 0, 0, 0)
#

def edge_walker(el, ose) :
    #print 'len(el) =', len(el)
    #print 'el=', el
    sp = None
    cnt = -1
    while( sp != ose ) :
        cnt += 1
        if( cnt == 0 ) :
            sp = ose
        pass
        fst = filter(lambda f: f[0:len(sp)] == sp, el)
        sp = fst[0][len(sp):]
        yield fst[0]
    pass
pass

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
##    ___  ____  _____  ________
##   / _ \/ __ \/  _/ |/ /_  __/
##  / ___/ /_/ // //    / / /
## /_/   \____/___/_/|_/ /_/


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
        self.tags = 0
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

    #----------------------------------------------------------------------------------------------------

    def calc_avg_parameter_data_from_point_list(self, gl) :
        #print '>>> TOP >>> ', here(), '  NAME = ',self.name
        
        #super(grid, self).__init__(0.0, 0.0, 0.0)
        self.__init__(0.0, 0.0, 0.0)

        ng = len(gl)
        if( ng == 0 ) : return

        
        #print 'POINT:', self

        #print 'GL = ', gl
        # AVERAGE X, Y, Z
        for g in gl :
            #print 'TYPE self = ', type(self), '   TYPE g =', type(g)
            for dd in range(len(self.v)) :
                self.v[dd] += g.v[dd]
            pass
        pass
        for dd in range(len(self.v)) :
             self.v[dd] /= float(ng) 
        pass

        # AVERAGE ANCILLARY DATA - IF POSSIBLE
        pd = {}

        # SUM
        for g in gl :
            #print 'LOOKING AT GRID PARAMS FOR GRID;', g
            for k, v in g.params.items() :
                #print '  LOOKING AT PARAM :', k
                if( k not in pd ) :
                    pd[k] = [v, 1]  # [ sum, count ] pairs
                    #print 'FIRST TIME FOR PARAM  VAL=', v
                else :
                    if( not isinstance(v, str) ) : # CAN'T SUM STRING PARAMS
                        pd[k][0] += v
                        pd[k][1] += 1
                        #print pd[k][1], ' TIME FOR PARAM  VAL=', v, '   RUNNING TOT = ', pd[k][0] 
                    pass
                pass
            pass
        pass

        # NOW AVERAGE
        for k, v  in pd.items() :
            if( not isinstance(v[0], str) ) :
                self.set_param(k, v[0] / float(v[1]))
            else :
               self.set_param(k, v[0])
            pass
            #print 'PARAM ', k, '  AVG VAL = ', self.get_param(k, -1.0)
        pass

    pass
        


    #----------------------------------------------------------------------------------------------------

pass



##########################################################################################################
# http://stackoverflow.com/questions/682504/what-is-a-clean-pythonic-way-to-have-multiple-constructors-in-python
##  _   _____________________  ___
## | | / / __/ ___/_  __/ __ \/ _ \
## | |/ / _// /__  / / / /_/ / , _/
## |___/___/\___/ /_/  \____/_/|_|
#
# ALSO USED AS AN EDGE OF A FACE  <-= NOTE
#
class vector(object) :

#class Zheese():
#    def __init__(self, *args, **kwargs):
#        #args -- tuple of anonymous arguments
#        #kwargs -- dictionary of named arguments
#        self.num_holes = kwargs.get('num_holes',random_holes())
        
#    def __init__ (self, x1, y1, z1, x2, y2, z2) :
    def __init__ (self, p1, p2 = None, parent = None) :
        self.tags = 0
        
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

        self.parent = parent

        self.calc_dv()

        #print 'dv =', self.dv
        #print 'LEN = ', self.len
        #print 'uv =', self.uv, np.linalg.norm(self.uv)
        
    pass

    #----------------------------------------------------------------------------------------------------

    def __repr__(self) :
        s = 'VECTOR: \n FROM: ' + str(self.sp) + '\n TO:   ' + str(self.ep) 
        s += '\n  DV = ' + str(self.dv) + '   LEN = ' + str(self.len) + '   UNIT_VEC = ' + str(self.uv) + '\n'
        return(s)
    pass
    
    #----------------------------------------------------------------------------------------------------

    def calc_dv(self) :
        self.dv = np.subtract(self.ep.v, self.sp.v)
        self.len = np.linalg.norm(self.dv)
        
        self.uv = np.array( [0.0, 0.0, 0.0] )
        if( self.len > 0.0 ) :
            self.uv = np.array( [ self.dv[0] / self.len,
                                  self.dv[1] / self.len,
                                  self.dv[2] / self.len ] )
        pass
    pass


    #----------------------------------------------------------------------------------------------------
    # FLIP IN PLACE
    def flip(self) :
        # SWAP START AND END POINTS
        self.sp, self.ep = self.ep, self.sp
        self.calc_dv()
    pass
    
    #----------------------------------------------------------------------------------------------------
    # FLIP A COPY
    def vector_from_flip(self) :
        v = vector(self.sp, self.ep, self.parent)
        v.flip()
        ## OR WE COULD HAVE JUST DID: v = vector(self.ep, self.sp)
        return(v)
    pass
    
    #----------------------------------------------------------------------------------------------------
    # JUST SWAP THE START AND END POINT LOCATIONS BUT LEAVE THE PARAMETERS THE SAME
    def vector_from_swap_geom_points(self) :
        v = vector(self.sp, self.ep, self.parent)
        v.sp.v, v.ep.v =  v.ep.v, v.sp.v
        v.calc_dv()
        return(v)
    pass
    
    #----------------------------------------------------------------------------------------------------
    
    def vec_cross(self, v) : 
        npw = np.cross( self.dv, v.dv ) 
        return(npw)
    pass

    #----------------------------------------------------------------------------------------------------

    def cross(self, v) : 
        npw = np.cross( self.dv, v.dv )
        nep = point(self.sp.v[0] + npw[0], self.sp.v[1] + npw[1], self.sp.v[2] + npw[2])
        for k, v  in self.sp.params.items() :
            nep.params[k] = v
        pass
        vv = vector(self.sp, nep, self.parent)
        return(vv)
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
    # SPECIFY A SCALE FACTOR - RETURNS A NEW VECTOR WITH THE SAME START POINT
    def vector_from_scale(self, sf) :
        rs = self.scale(sf)
        rest = OrderedDict( [(k, rs[k]) for k in rs.keys()[3:] ] )
        p = point(rs['X'], rs['Y'], rs['Z'], rest)
        vv = vector(self.sp, p, self.parent)
        return(vv)
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
    #----------------------------------------------------------------------------------------------------

    # SPECIFY A LENGTH - RETURNS A NEW VECTOR WITH SAME START POINT - DEFAULTS TO A UNIT VECTOR
    def vector_from_length(self, new_len = 1.0) :
        sf = 1.0
        if( self.len != 0.0 ) :
            sf = new_len / self.len
        pass
        p =  self.point_from_scale(sf)
        vv = vector(self.sp, p, self.parent)
        return(vv)
    pass

pass

##########################################################################################################
##    _______  _________
##   / __/ _ |/ ___/ __/
##  / _// __ / /__/ _/
## /_/ /_/ |_\___/___/


class face(object) :

    #----------------------------------------------------------------------------------------------------

    def __init__ (self, name) :
        self.tags = 0
        self.name = name
        self.edges = []
        #print "nedges passed = ", len(el) 
    pass

    #----------------------------------------------------------------------------------------------------

    def define_face_from_edge_list(self, *el) :
        for i, edg in enumerate(el) :
            #print 'ADDING FACE EDGE # ', i
            self.edges.append(copy.deepcopy(edg)) 
        pass
        for e in self.edges :
            e.parent = self
        pass
    pass

    #----------------------------------------------------------------------------------------------------
    # NORMALS POINT IN
    def define_face_from_points(self, p1, p2, p3, p4) :
        self.edges = []
        self.edges.append( vector(p1, p2, self) )
        self.edges.append( vector(p2, p3, self) )
        self.edges.append( vector(p3, p4, self) )
        self.edges.append( vector(p4, p1, self) )
    pass

    #----------------------------------------------------------------------------------------------------

    def get_face_corner_points(self) :
        cps = []
        cps.append(self.edges[0].sp)
        cps.append(self.edges[1].sp)
        cps.append(self.edges[2].sp)
        cps.append(self.edges[3].sp)
        return(cps)
    pass

    #----------------------------------------------------------------------------------------------------

    def __del__ (self) :
        del self.edges
    pass

    #----------------------------------------------------------------------------------------------------

    def __repr__(self):
        s = 'FACE:  ' + self.name + '  EDGES...\n'
        #print 'NUMBER OF EDGES =', len(self.edges)
        for i, e in enumerate(self.edges) :
            s += '>>> ' + str(i) + ') ' + str(e) + ' <<< ps=' + str(e.sp.get_param('PRES', 0.0)) + ' pe=' + str(e.ep.get_param('PRES', 0.0)) + ' \n'
        pass
        return(s)
    pass

    #----------------------------------------------------------------------------------------------------

pass

##########################################################################################################
##    ___  __   ____  _______ __
##   / _ )/ /  / __ \/ ___/ //_/
##  / _  / /__/ /_/ / /__/ ,<
## /____/____/\____/\___/_/|_|

# SPECIAL CASE NEEDED - ONE FACE THAT REPRESENTS THE MIDCAMBER
#   THEN CREATE THE 3D BLOCK FROM THAT MIDCAMBER SURFACE (WITH THICKNESSES)

# BLOCK :  +/- 1 parametric CUBE CENTERED AT ZERO
# EDGE INDICIES (su, sv, sw, eu, ev, ew) ALWAYS POINT FROM LO TO HI VALUE
# FIRST THREE INDICIES ARE THE START POINT THE LAST THREE ARE THE END POINT INDICIES
#  e.g.   EDGE(-1, 1, 1, 1, 1, 1) POINTS FROM POINT (-1, 1, 1) TO (1, 1, 1)
#     THE EDGE VECTOR POINTS DOWN THE U AXIS ALONG THE EDGE WHERE V, ARE W ARE THIER MAXIMUM
#  ONLY ONE OF THE INDICIES WILL CHANGE FROM (-) TO (+) FROM THE START POINT TO THE END POINT
#  EACH BLOCK STORES 8 CORNER POINTS AND 12 EDGE VECTORS
    
##                            W
##                          |                   .-' V
##      -1,1,1              |                .-'
##               .:.........|..............-'1,1,1
##             .' |         |           .'
## -1,-1,1   .:....................._.-' |
##           |    |         |     .-| 1,-1,1
##           |    |         |  .-'  |    |
##           |    |         .-:.....|.......................U
##           |    |      0,0,0      |    |
##           |    |                 |    |
##           |    |-1,1,-1          |    |
##           |   .'.................|....'   1, 1, -1
##           |  /                   | .'
##           `.:....................|'
##
##        -1,-1,-1                 1,-1,-1

class block(object) :



    #----------------------------------------------------------------------------------------------------

    def __init__ (self, name, model) :
        self.tags = 0
        self.defined_by_tag = -1
        self.name = name

        self.mesh = None

        self.model = model
        self.faces = {} # FACES DICT
        self.pl = {} # POINT DICT
        self.el = {} # EDGE DICT
        
        self.min_ne = 2 # min number of elements in any direction
        self.avg_dv = [ np.zeros(3) ] * 3
        self.coords = [ np.zeros(3) ] * 3
        self.approx_middle  = np.zeros(3)
        self.edge_sizes  = np.zeros(3) # MIN, AVG, MAX
        self.geom_jar = np.zeros((3, 3)) # TBD
        self.nom_el_size = np.zeros(3)
        self.ne = np.zeros(3, dtype=np.int32)
        self.ng = np.zeros(3, dtype=np.int32)
        # MIGHT WANT TO TWEEK THIS (WITHING ACCEPTABLE VALUES) TO ACCOMMODATE DISSIMILAR MESHES
        # 1/40 < TARGET_AR < 40  ???  
        self.target_ar = None
        self.fne = np.zeros(3, dtype=np.int32)


    pass

    #----------------------------------------------------------------------------------------------------

    def __del__ (self) :
        del self.faces
        del self.mesh
        del self.pl
        del self.el
    pass

    #----------------------------------------------------------------------------------------------------

    def __repr__(self) :
        s = 'BLOCK:  FACES...\n'
        print 'NUMBER OF FACES =', len(self.faces)
        for i, f in enumerate(self.faces.values()) :
            s += '>>> ' + str(i) + ') ' + str(f.name) + ' <<<\n'
        pass
        
        
        return(s)
    pass

    #----------------------------------------------------------------------------------------------------

    def force_element_count(self, tne) :
        if( isinstance(tne, list) ) :
            i = -1
            for n in tne :
                i += 1
                if( n is None ) : continue
                if( n > 0 ) :
                    self.fne[i] = n
                pass
            pass
        else :
            for i in range(3) :
                self.fne[i] = tne
            pass
        pass
    pass
        

    #----------------------------------------------------------------------------------------------------

    def get_grid_list_from_tags(self, tag, igl = None) :
        rgl = self.mesh.get_grid_list_from_tags(tag, igl)
        return(rgl)
    pass

    #----------------------------------------------------------------------------------------------------

    def get_grid_list_from_not_tags(self, tag, igl = None) :
        rgl = self.mesh.get_grid_list_from_not_tags(tag, igl)
        return(rgl)
    pass

    #----------------------------------------------------------------------------------------------------

        
    def set_material(self, mat) :
        self.mesh.set_material(mat)
    pass
    
    #----------------------------------------------------------------------------------------------------

    # FROM CALCULIX MANUAL...
    #  face 1: 1-2-3-4 -> FTAGS.MIN_W - 1
    #  face 2: 5-8-7-6 -> FTAGS.MAX_W - 2
    #  face 3: 1-5-6-2 -> FTAGS.MIN_V - 4
    #  face 4: 2-6-7-3 -> FTAGS.MAX_U - 8
    #  face 5: 3-7-8-4 -> FTAGS.MAX_V - 16
    #  face 6: 4-8-5-1 -> FTAGS.MIN_U - 32
    #
    # FACE NORMALS POINT IN
    #
##                            
##                            W
##                          |                   .-' V
##      -1,1,1              |                .-'
##               .:.........|..............-'1,1,1
##             .' |         |           .'
## -1,-1,1   .:....................._.-' |
##           |    |         |     .-| 1,-1,1
##           |    |         |  .-'  |    |
##           |    |         .-:.....|.......................U
##           |    |      0,0,0      |    |
##           |    |                 |    |
##           |    |-1,1,-1          |    |
##           |   .'.................|....'   1, 1, -1
##           |  /                   | .'
##           `.:....................|'
##
##        -1,-1,-1                 1,-1,-1
##                                                   
##                8        +W               7    
##               .:........................  
##             .' |      F2             .'
##         5 .:....................._.-' |
##           |    |                 |6   |
##   -U      |    |      F5(back)   |    |      +U
##           |    |                 |    |                   
##           | F6 |    F3(front)    | F4 |
##           |    |                 |    |
##           |    |4                |    |
##           |   .'.................|....' 3
##           |  /        F1         | .'
##           `.:....................|'
##           1                      2
##                     -W                            

    #----------------------------------------------------------------------------------------------------
    # General block definiton by points
    #    - We use points here therefore we can NOT change pressures independently on each face  

    def define_block_from_points(self, p1, p2, p3, p4, p5, p6, p7, p8) :
        self.faces.clear()
        
        self.faces[FTAGS.MIN_W] = face('F1')
        self.faces[FTAGS.MIN_W].define_face_from_points(p1, p2, p3, p4)

        self.faces[FTAGS.MAX_W] = face('F2')
        self.faces[FTAGS.MAX_W].define_face_from_points(p5, p8, p7, p6)

        self.faces[FTAGS.MIN_V] = face('F3')
        self.faces[FTAGS.MIN_V].define_face_from_points(p1, p5, p6, p2)
        
        self.faces[FTAGS.MAX_U] = face('F4')
        self.faces[FTAGS.MAX_U].define_face_from_points(p2, p6, p7, p3)
        
        self.faces[FTAGS.MAX_V] = face('F5')
        self.faces[FTAGS.MAX_V].define_face_from_points(p3, p7, p8, p4)
        
        self.faces[FTAGS.MIN_U] = face('F6')
        self.faces[FTAGS.MIN_U].define_face_from_points(p4, p8, p5, p1)


        self.defined_by_tag = BTAGS.POINTS
        self.add_points_and_edges_from_faces()

    pass

    #----------------------------------------------------------------------------------------------------
    # General block definiton by faces
    #    - We use faces here (not points) so we CAN change pressures independently on each face  

    def define_block_from_faces(self, f1, f2, f3, f4, f5, f6) :
        self.faces.clear()
        
        # WE NEED TO STORE THE CORRESPONDING POINTS AND EDGES SO DO A DEEP COPY
        self.faces[FTAGS.MIN_W] = copy.deepcopy(f1)
        self.faces[FTAGS.MAX_W] = copy.deepcopy(f2)
        self.faces[FTAGS.MIN_V] = copy.deepcopy(f3)
        self.faces[FTAGS.MAX_U] = copy.deepcopy(f4)
        self.faces[FTAGS.MAX_V] = copy.deepcopy(f5)
        self.faces[FTAGS.MIN_U] = copy.deepcopy(f6)

        self.defined_by_tag = BTAGS.FACES
        self.add_points_and_edges_from_faces()

    pass

    #----------------------------------------------------------------------------------------------------
    # RIGHT HAND RULE IN POINT LISTING
    def define_block_from_mid_plane_points(self, p1, p2, p3, p4) :
        self.faces.clear()       
        ff = face(self.name + '_MIDSURF')
        ff.define_face_from_points(p1, p2, p3, p4)
        ff.tags |= FTAGS.MID_W
        self.defined_by_tag = BTAGS.EXPLODED
        self.explode_block_from_middle_surface(ff)
    pass

    #----------------------------------------------------------------------------------------------------

    def explode_block_from_middle_surface(self, f) :

        # CURRENT EDGE INDEX
        for cei in (3, 0, 1, 2) :
            nei = cei + 1
            
            if( nei == len(f.edges) ) :
                nei = 0
            pass
        
            #print 'CEI =', cei, '  NEI =', nei

            # THIS IS THE COMMON POINT U,V INDEX
            # THE FOLLOWING ASSUMES 4 EDGES PER FACE
            u = -1
            if( ( cei == 0 ) or ( cei == 1 ) ) : u = 1
            v = -1
            if( ( cei == 1 ) or ( cei == 2 ) ) : v = 1


            nu = -1
            if( ( nei == 0 ) or ( nei == 1 ) ) : nu = 1
            nv = -1
            if( ( nei == 1 ) or ( nei == 2 ) ) : nv = 1

            #print 'LEN EDGES = ', len(f.edges)
            ce = f.edges[cei]
            #print 'CURRENT EDGE=', ce
            w = 0
            self.pl[u, v, w] = ce.ep

            #print 'ADDED POINT MID u=', u, '  v=', v, '  w=', w, '   ', self.pl[u,v,w]
            
            ne = f.edges[nei]

            #print u, v, 0, '-=>', nu, nv, 0

            #self.el[u, v, 0, nu, nv, 0] = ne
  
            # FLIP THE CURRENT VECTOR
            ce_flipped = ce.vector_from_swap_geom_points()
            #print 'ce_flipped = ', ce_flipped
            #ce_flipped = ce.vector_from_scale(-1.0)
            # GET THE NORMAL AND SCALE IT TO MATCH THE HALF THICKNESS
            th = ne.sp.get_param('THICK', 1.0)
            nv = ne.cross(ce_flipped).vector_from_length(th/2.0)
            w = 1
            self.pl[u, v, w] = nv.ep
            #print 'ADDED POINT TOP u=', u, '  v=', v, '  w=', w, '   ', self.pl[u,v,w]

            # FLIP THE NORMAL TO GET OTHER CORNER POINT
            # WITH vector_from_scale(-1.0) THE START POINT STAYS THE SAME
            nnv = nv.vector_from_scale(-1.0)
            #print 'NV=', nv
            #print 'NNV=', nnv
            w = -1
            self.pl[u, v, w] = nnv.ep
            #print 'ADDED POINT BOT u=', u, '  v=', v, '  w=', w, '   ', self.pl[u,v,w]

            #tedge = vector(self.pl[u, v, -1], self.pl[u, v, 1], self )
            #self.el[u, v, -1, u, v, 1] = tedge
            
        pass

        self.add_edges_and_faces_from_points()

        #print '_'*80 + '\n'
        self.print_points()
        #print '_'*80 + '\n'
        
        self.print_edges()

    pass

    #----------------------------------------------------------------------------------------------------
    # THIS IS ONLY CALLED FOR FACE DEFINED BLOCKS
    # IN THE GENERAL CASE WE CAN HAVE MORE THAN ONE POINT AT A CORNER
    #  THIS IS SO WE CAN DE DEFINE NEIGHBORING FACES WITH DIFFERENT POINTS WITH DIFFERENT PRESSURES, ETC.
    # WE MUST ASSUME THAT THE EDGES AND THE POINTS DEFINED BY THE FACES ARE COINCIDENT OTHERWISE
    #  WE CAN'T MESH THE BLOCK
    def add_points_and_edges_from_faces(self) :
        
        p = 1
        n = -1

        # FACE NORMALS POINT IN
        
        f = self.faces[FTAGS.MIN_W]

        # nnn is defined by
        #  FACE F1, F3, F6
        
        
        self.pl[n, n, n] = copy.deepcopy(f.edges[0].sp)
        self.pl[p, n, n] = copy.deepcopy(f.edges[0].ep)
        self.pl[p, p, n] = copy.deepcopy(f.edges[2].sp)
        self.pl[n, p, n] = copy.deepcopy(f.edges[2].ep)
        
        self.el[n, n, n, p, n, n] = copy.deepcopy(f.edges[0])
        self.el[p, n, n, p, p, n] = copy.deepcopy(f.edges[1])
        self.el[n, p, n, p, p, n] = copy.deepcopy(f.edges[2].vector_from_flip())
        self.el[n, n, n, n, p, n] = copy.deepcopy(f.edges[3].vector_from_flip())


        f = self.faces[FTAGS.MAX_W]
        self.pl[n, n, p] = copy.deepcopy(f.edges[0].sp)
        self.pl[n, p, p] = copy.deepcopy(f.edges[0].ep)
        self.pl[p, p, p] = copy.deepcopy(f.edges[2].sp)
        self.pl[p, n, p] = copy.deepcopy(f.edges[2].ep)

        self.el[n, n, p, n, p, p] = copy.deepcopy(f.edges[0])
        self.el[n, p, p, p, p, p] = copy.deepcopy(f.edges[1])
        self.el[p, n, p, p, p, p] = copy.deepcopy(f.edges[2].vector_from_flip())
        self.el[n, n, p, p, n, p] = copy.deepcopy(f.edges[3].vector_from_flip())


        f = self.faces[FTAGS.MIN_V]
        self.el[n, n, n, n, n, p] = copy.deepcopy(f.edges[0])
        self.el[p, n, n, p, n, p] = copy.deepcopy(f.edges[2].vector_from_flip())
        
        f = self.faces[FTAGS.MAX_V]
        self.el[p, p, n, p, p, p] = copy.deepcopy(f.edges[0])
        self.el[n, p, n, n, p, p] = copy.deepcopy(f.edges[2].vector_from_flip())


        self.calc_geom_jar()
        
        pass
    pass

    #----------------------------------------------------------------------------------------------------
    # THIS IS ONLY CALLED FOR EXLODED BLOCKS
    def add_edges_and_faces_from_points(self) :

        p = 1
        n = -1

        # NOTE:
        # The edges stored in the block class "el" dictionary always point
        # from the minimum parameter direction to the maximum parameter direction
        #
        # The edges stored in the face are stored so they are oriented so they form a loop
        # that circumventes the perimeter of the face, eventually ending back to where it started.
        
        #self.avg_dv = [ np.zeros(3) ] * 3
        #self.coords = [ np.zeros(3) ] * 3
        #half_dv = np.zeros(3)
        
        #print '\nAVG_DV[]=', self.avg_dv
       
        #c = 0
        # EDGES THAT GO FROM LOU TO HIU
        for si in [ (n, n, n), (n, p, n), (n, p, p), (n, n, p) ] :
            #c += 1
            # CHANGE ONE INDEX IN THE DIRECTION OF CONCERN TO GET END POINT INDEX
            lei = list(si) ; lei[DIR.U] = p ; ei = tuple(lei)
            edge_index = si + ei        
            self.el[edge_index] = vector(self.pl[si], self.pl[ei], self)
            #print 'DELU ', si, ' -=> ', ei, ' DV=', self.el[edge_index].dv   
            #self.avg_dv[DIR.U] = self.avg_dv[DIR.U] + self.el[edge_index].dv
        pass
        #print DIR.U, 'DELU AVG = ',  self.avg_dv[DIR.U]   
        #self.avg_dv[DIR.U] = self.avg_dv[DIR.U] / float(c)
        #print DIR.U, 'DELU AVG = ',  self.avg_dv[DIR.U]   

        #print '\nAVG_DV[]=', self.avg_dv
        
        #c = 0
        #print 'PRE:', DIR.V, 'DELV AVG = ',  self.avg_dv[DIR.V]   
        # EDGES THAT GO FROM LOV TO HIV
        for si in [ (n, n, n), (p, n, n), (p, n, p), (n, n, p) ] :
            #c += 1
            # CHANGE ONE INDEX IN THE DIRECTION OF CONCERN TO GET END POINT INDEX
            lei = list(si) ; lei[DIR.V] = p ; ei = tuple(lei)
            edge_index = si + ei
            #print 'PSTART ', self.pl[si]
            #print 'PEND ', self.pl[ei]
            self.el[edge_index] = vector(self.pl[si], self.pl[ei], self)
            #print 'BEFORE ADD =', self.avg_dv[DIR.V]
            #print 'DELV ', si, ' -=> ', ei, ' DV=', self.el[edge_index].dv   
            #self.avg_dv[DIR.V] = self.avg_dv[DIR.V] +self.el[edge_index].dv
            #print 'AFTER ADD =', self.avg_dv[DIR.V]
        pass
        #print DIR.V, 'DELV AVG = ',  self.avg_dv[DIR.V]   
        #self.avg_dv[DIR.V] = self.avg_dv[DIR.V] / float(c)
        #print DIR.V, 'DELV AVG = ',  self.avg_dv[DIR.V]   

        #print '\nAVG_DV[]=', self.avg_dv

        #c = 0
        # EDGES THAT GO FROM LOW TO HIW
        for si in [ (n, n, n), (p, n, n), (p, p, n), (n, p, n) ] :
            #c += 1
            # CHANGE ONE INDEX IN THE DIRECTION OF CONCERN TO GET END POINT INDEX
            lei = list(si) ; lei[DIR.W] = p ; ei = tuple(lei)
            edge_index = si + ei        
            self.el[edge_index] = vector(self.pl[si], self.pl[ei], self)
            #print 'DELW ', si, ' -=> ', ei, ' DV=', self.el[edge_index].dv   
            #self.avg_dv[DIR.W] = self.avg_dv[DIR.W] + self.el[edge_index].dv
        pass
        #print DIR.W, 'DELV AVG = ',  self.avg_dv[DIR.V]   
        #self.avg_dv[DIR.W] = self.avg_dv[DIR.W] / float(c)
        #print DIR.W, 'DELV AVG = ',  self.avg_dv[DIR.V]   

        # NORMALS POINT IN
        #  face 1: 1-2-3-4 -> FTAGS.MIN_W - 1
        #  face 2: 5-8-7-6 -> FTAGS.MAX_W - 2
        #  face 3: 1-5-6-2 -> FTAGS.MIN_V - 4
        #  face 4: 2-6-7-3 -> FTAGS.MAX_U - 8
        #  face 5: 3-7-8-4 -> FTAGS.MAX_V - 16
        #  face 6: 4-8-5-1 -> FTAGS.MIN_U - 32
        #
        #pordertoindex = [ (n, n, n), (p, n, n), (p, p, n), (n, p, n), (n, n, p), (p, n, p), (p, p, p), (n, p, p) ]

        # FACE 1 - MIN_W
        edges = []
        edges.append(self.el[n, n, n, p, n, n])
        edges.append(self.el[p, n, n, p, p, n])
        edges.append(self.el[n, p, n, p, p, n].vector_from_flip())
        edges.append(self.el[n, n, n, n, p, n].vector_from_flip())
        ff = face(self.name + '_MIN_W', *edges)
        ff.tags |= FTAGS.MIN_W
        self.faces[FTAGS.MIN_W] = ff

        # FACE 2 - MAX_W
        edges = []
        edges.append(self.el[n, n, p, p, n, p])
        edges.append(self.el[p, n, p, p, p, p])
        edges.append(self.el[n, p, p, p, p, p].vector_from_flip())
        edges.append(self.el[n, n, p, n, p, p].vector_from_flip())
        ff = face(self.name + '_MAX_W', *edges)
        ff.tags |= FTAGS.MAX_W
        self.faces[FTAGS.MAX_W] = ff
        
        # FACE 3 - MIN_V
        edges = []
        edges.append(self.el[n, n, n, p, n, n])
        edges.append(self.el[p, n, n, p, n, p])
        edges.append(self.el[n, n, p, p, n, p].vector_from_flip())
        edges.append(self.el[n, n, n, n, n, p].vector_from_flip())
        ff = face(self.name + '_MIN_V', *edges)
        ff.tags |= FTAGS.MIN_V
        self.faces[FTAGS.MIN_V] = ff

        # FACE 4 - MAX_U
        edges = []
        edges.append(self.el[p, n, n, p, p, n])
        edges.append(self.el[p, p, n, p, p, p])
        edges.append(self.el[p, n, p, p, p, p].vector_from_flip())
        edges.append(self.el[p, n, n, p, n, p].vector_from_flip())
        ff = face(self.name + '_MAX_U', *edges)
        ff.tags |= FTAGS.MAX_U
        self.faces[FTAGS.MAX_U] = ff

        # FACE 5 - MAX_V
        edges = []
        edges.append(self.el[n, p, n, p, p, n])
        edges.append(self.el[p, p, n, p, p, p])
        edges.append(self.el[n, p, p, p, p, p].vector_from_flip())
        edges.append(self.el[n, p, n, n, p, p].vector_from_flip())
        ff = face(self.name + '_MAX_V', *edges)
        ff.tags |= FTAGS.MAX_V
        self.faces[FTAGS.MAX_V] = ff

        # FACE 6 - MIN_U
        edges = []
        edges.append(self.el[n, n, n, n, p, n])
        edges.append(self.el[n, p, n, n, p, p])
        edges.append(self.el[n, n, p, n, p, p].vector_from_flip())
        edges.append(self.el[n, n, n, n, n, p].vector_from_flip())
        ff = face(self.name + '_MIN_U', *edges)
        ff.tags |= FTAGS.MIN_U
        self.faces[FTAGS.MIN_U] = ff


        ## print 'FACES FROM EXPLOSION...'
        ## for f in self.faces.values() :
        ##     print f
        ## pass

        #print self.el

        #half_dv = ( self.avg_dv[DIR.U] / 2.0 ) + ( self.avg_dv[DIR.V] / 2.0 ) + ( self.avg_dv[DIR.W] / 2.0 )

        #print 'BLOCK BASE = ', self.pl[n, n, n].v
        #print 'AVG_DV = \n', self.avg_dv
        #print 'HALF_DV = ', half_dv
        
        # self.approx_middle = self.pl[n, n, n].v + half_dv
        #self.approx_middle = np.zeros(3)
        ## self.edge_sizes = [ constants.BIG_REAL, 0.0, -1.0 * constants.BIG_REAL ]  # MIN, AVG, MAX
        ## for kk, ev in self.el.items() :
        ##     self.edge_sizes[constants.MIN] = min(self.edge_sizes[constants.MIN], ev.len)
        ##     self.edge_sizes[constants.AVG] = self.edge_sizes[constants.AVG] + ev.len
        ##     self.edge_sizes[constants.MAX] = max(self.edge_sizes[constants.MAX], ev.len)
        ## pass
        ## self.edge_sizes[constants.AVG] = self.edge_sizes[constants.AVG] / float(len(self.el))

        ## #print 'MIN EDGE = ', self.edge_sizes[constants.MIN] 
        ## #print 'AVG EDGE = ', self.edge_sizes[constants.AVG] 
        ## #print 'MAX EDGE = ', self.edge_sizes[constants.MAX] 

        ## for kk, pv in self.pl.items() :
        ##     self.approx_middle =  self.approx_middle + pv.v
        ## pass
        ## self.approx_middle =  self.approx_middle / len(self.pl)
        ## #print 'APPROX MIDDLE = ', self.approx_middle

        ## # CALCULATE AN ORTHOGONAL BLOCK COODINATE SYSTEM
        ## self.coords[DIR.U] = self.avg_dv[DIR.U]
        ## self.coords[DIR.W] = np.cross( self.coords[DIR.U], self.avg_dv[DIR.V] )
        ## self.coords[DIR.V] = np.cross( self.coords[DIR.W], self.coords[DIR.U] )

        ## # NORMALIZE
        ## for r in range(3) :
        ##     self.coords[r] = self.coords[r] / np.linalg.norm(self.coords[r])
        ## pass

        
        self.calc_geom_jar()
        
    pass

    #----------------------------------------------------------------------------------------------------
    # BLOCK GEOMETRY JACOBIAN ASPECT RATIO (JAR)
    def calc_geom_jar(self) :
        p = 1
        n = -1


        ## 'BLOCK EDGE KEYS...'
        ## for k, v in self.el.items() :
        ##     print k
        ## pass




        # AVG ALL THE FACE VALUES TO GET NEW BLOCK.EDGE AND BLOCK.POINT VALUES
        # P1 - start/end pairs OF FACE EDGES
        til = [ [ 1, (n,n,n), [ [ FTAGS.MIN_U, [ 3, 2 ] ],
                                [ FTAGS.MIN_V, [ 0, 3 ] ],
                                [ FTAGS.MIN_W, [ 0, 3 ] ],
                                ] ],
                [ 2, (p,n,n), [ [ FTAGS.MAX_U, [ 0, 3 ] ],
                                [ FTAGS.MIN_V, [ 3, 2 ] ],
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

        kl = self.el.keys()
        tp = point()
        pl = []

        for pn, pi, fiv in til :
            #print 'pn=', pn
            #print 'pi=', pi
            #print 'fiv=', fiv
            tp.__init__()
            del pl[:]
            for ft, il in fiv:
                pl.append(self.faces[ft].edges[il[0]].sp)
                pl.append(self.faces[ft].edges[il[1]].ep)
                #print 'FACE TAG =', ft,
                #print '  START POINT EDGE INDEX =', il[0], '  PRESSURE =', self.faces[ft].edges[il[0]].sp.get_param('PRES', 0.0)
                #print '  END POINT EDGE INDEX =', il[1], '  PRESSURE =', self.faces[ft].edges[il[1]].sp.get_param('PRES', 0.0)
            pass
            tp.calc_avg_parameter_data_from_point_list(pl)
            self.pl[pi].v = copy.deepcopy(tp.v)
            self.pl[pi].params = copy.deepcopy(tp.params)

            
            # NOW WE NEED TO MAKE SURE ALL BLOCK EDGES REFERNCE THE BLOCK CORNER POINT
            fkl_start = filter( ( lambda k: k[0:3] == pi ), kl)
            for k in fkl_start :
                self.el[k].sp = copy.deepcopy(self.pl[pi])
                self.el[k].calc_dv()
            pass
            fkl_end = filter( ( lambda k: k[3:] == pi ), kl)
            for k in fkl_end :
                self.el[k].ep = copy.deepcopy(self.pl[pi])
                self.el[k].calc_dv()
            pass

        pass


        # CALCULATE NOMINAL BLOCK DIMENSIONS AND COORDIANTES
        c = 0
        # EDGES THAT GO FROM LOU TO HIU
        for si in [ (n, n, n), (n, p, n), (n, p, p), (n, n, p) ] :
            c += 1
            # CHANGE ONE INDEX IN THE DIRECTION OF CONCERN TO GET END POINT INDEX
            lei = list(si) ; lei[DIR.U] = p ; ei = tuple(lei)
            edge_index = si + ei        
            self.avg_dv[DIR.U] = self.avg_dv[DIR.U] + self.el[edge_index].dv
        pass
        self.avg_dv[DIR.U] = self.avg_dv[DIR.U] / float(c)

        c = 0
        # EDGES THAT GO FROM LOV TO HIV
        for si in [ (n, n, n), (p, n, n), (p, n, p), (n, n, p) ] :
            c += 1
            # CHANGE ONE INDEX IN THE DIRECTION OF CONCERN TO GET END POINT INDEX
            lei = list(si) ; lei[DIR.V] = p ; ei = tuple(lei)
            edge_index = si + ei
            self.avg_dv[DIR.V] = self.avg_dv[DIR.V] + self.el[edge_index].dv
        pass
        self.avg_dv[DIR.V] = self.avg_dv[DIR.V] / float(c)

        c = 0
        # EDGES THAT GO FROM LOW TO HIW
        for si in [ (n, n, n), (p, n, n), (p, p, n), (n, p, n) ] :
            c += 1
            # CHANGE ONE INDEX IN THE DIRECTION OF CONCERN TO GET END POINT INDEX
            lei = list(si) ; lei[DIR.W] = p ; ei = tuple(lei)
            edge_index = si + ei        
            self.avg_dv[DIR.W] = self.avg_dv[DIR.W] + self.el[edge_index].dv
        pass
        self.avg_dv[DIR.W] = self.avg_dv[DIR.W] / float(c)


        self.edge_sizes = [ constants.BIG_REAL, 0.0, -1.0 * constants.BIG_REAL ]  # MIN, AVG, MAX
        for kk, ev in self.el.items() :
            self.edge_sizes[constants.MIN] = min(self.edge_sizes[constants.MIN], ev.len)
            self.edge_sizes[constants.AVG] = self.edge_sizes[constants.AVG] + ev.len
            self.edge_sizes[constants.MAX] = max(self.edge_sizes[constants.MAX], ev.len)
        pass
        self.edge_sizes[constants.AVG] = self.edge_sizes[constants.AVG] / float(len(self.el))

        #print 'MIN EDGE = ', self.edge_sizes[constants.MIN] 
        #print 'AVG EDGE = ', self.edge_sizes[constants.AVG] 
        #print 'MAX EDGE = ', self.edge_sizes[constants.MAX] 

        for kk, pv in self.pl.items() :
            self.approx_middle =  self.approx_middle + pv.v
        pass
        self.approx_middle =  self.approx_middle / len(self.pl)
        #print 'APPROX MIDDLE = ', self.approx_middle

        # CALCULATE AN ORTHOGONAL BLOCK COODINATE SYSTEM
        self.coords[DIR.U] = self.avg_dv[DIR.U]
        self.coords[DIR.W] = np.cross( self.coords[DIR.U], self.avg_dv[DIR.V] )
        self.coords[DIR.V] = np.cross( self.coords[DIR.W], self.coords[DIR.U] )

        # NORMALIZE
        for r in range(3) :
            self.coords[r] = self.coords[r] / np.linalg.norm(self.coords[r])
        pass


        # MIGHT WANT TO TWEEK THIS (WITHING ACCEPTABLE VALUES) TO ACCOMMODATE DISSIMILAR MESHES
        # 1/40 < TARGET_AR < 40  ???  
        if( self.target_ar is None ) :
            self.target_ar = np.ones(3)
        pass

        self.geom_jar = np.zeros((3, 3))
        self.nom_el_size = np.zeros(3)
        self.ne = np.zeros(3, dtype=np.int32)
        self.ne.fill(self.min_ne)
        self.ng = np.zeros(3, dtype=np.int32)


        dvs = np.zeros(3)
        dvs[DIR.U] =  np.linalg.norm(self.avg_dv[DIR.U])
        dvs[DIR.V] =  np.linalg.norm(self.avg_dv[DIR.V])
        dvs[DIR.W] =  np.linalg.norm(self.avg_dv[DIR.W])

        #print 'DVS = ', dvs

        for r in range(3) :
            for c in range(3) :
                self.geom_jar[r, c] = dvs[r] / dvs[c]
            pass
        pass

        #print 'DV_U = ', self.avg_dv[DIR.U], '  DV = ', dvs[DIR.U]
        #print 'DV_V = ', self.avg_dv[DIR.V], '  DV = ', dvs[DIR.V]
        #print 'DV_W = ', self.avg_dv[DIR.W], '  DV = ', dvs[DIR.W]
        
        #print 'GEOM JAR =\n', self.geom_jar

        
        id_short = dvs.argmin()
        #print 'id_short = ', id_short

        self.nom_el_size[id_short] = dvs[id_short] / self.ne[id_short]

        for i in range(3) :
            if( i == id_short ) : continue
            self.nom_el_size[i] = self.target_ar[i] * self.nom_el_size[id_short]
            self.ne[i] = dvs[i] / self.nom_el_size[i]
            if( (self.ne[i] % 2) == 1 ) :
                self.ne[i] = self.ne[i] - 1
                if( self.ne[i] < self.min_ne ) :
                    self.ne[i] = self.min_ne
                pass
                self.nom_el_size[i] = dvs[i] / self.ne[i]
            pass
        pass


        ## # DT test
        ## self.ne[0] = 4
        ## self.ne[1] = 4
        ## self.ne[2] = 2

        ## if( self.name == 'A' ) :
        ##     self.ne[0] = 4
        ##     self.ne[1] = 4
        ##     self.ne[2] = 4
        ## pass


        # USER OVERIDE ELEMENT COUNT
        print 'FNE =', self.fne
        print 'NE =', self.ne
        
        for i in range(3) :
            if( self.fne[i] > 0 ) :
                self.ne[i] = self.fne[i]
            pass
        pass
        print 'NE =', self.ne
        

        kl = self.faces.keys()
        kl.sort()
        for k in kl :
            print self.faces[k]
        pass



        print 'NE = \n', self.ne 
        print 'NOM_EL_SIZE = \n', self.nom_el_size 

    pass

    #----------------------------------------------------------------------------------------------------

    def print_points(self) :
        
        kl = sorted(self.pl.keys())
        print kl
        
        print '_'*80 + '\n'
        print 'POINTS: BLOCK: |' + self.name + '| # OF POINTS = ', len(self.pl)
        for k in kl :
            print k, ' -=>  ', self.pl[k]
        pass
        print '_'*80 + '\n'
        sys.stdout.flush()
    pass

    #----------------------------------------------------------------------------------------------------

    def print_edges(self) :
        
        kl = sorted(self.el.keys())

        print '_'*80 + '\n'
        print 'EDGES: BLOCK: |' + self.name + '| # OF EDGES = ', len(self.el)
        for k in kl :
            print k, ' -=>  ', self.el[k]
        pass
        print '_'*80 + '\n'
        sys.stdout.flush()
    pass

    #----------------------------------------------------------------------------------------------------

    def plot_block_on_ax(self, ax) :
        
        sys.stdout.flush()

        # USED EDGE LIST TO PRINT THE BLOCKS
        # GET EDGES WITH NO ZERO INDICES
        fel = filter(lambda f: np.prod(f) , self.el)
        #print 'FEL=', fel
        for pp in fel :
            
            
            pl = list( pp )
            #print 'pl=', pl
            xi = pl[0::3]
            yi = pl[1::3]
            zi = pl[2::3]

            #print 'XiYiZi=', xi, yi, zi

            
            p1 = self.pl[xi[0], yi[0], zi[0]]
            p2 = self.pl[xi[1], yi[1], zi[1]]

            x = [ p1.v[0], p2.v[0] ]
            y = [ p1.v[1], p2.v[1] ]
            z = [ p1.v[2], p2.v[2] ]
            
            ax.plot(x, y, z, zdir='z', color='black')
            #ax.plot(x, y, z)
        pass

        ## ax.text(self.approx_middle[0],
        ##         self.approx_middle[1],
        ##         self.approx_middle[2], 'O', color='red')

        ax.scatter(self.approx_middle[0],
                self.approx_middle[1],
                self.approx_middle[2], c='red', marker='o')

        ep = [ [ self.coords[DIR.U], 'RED' ],
               [ self.coords[DIR.V], 'GREEN' ],
               [ self.coords[DIR.W], 'BLUE'] ]

        ff = self.edge_sizes[constants.AVG] / 4.0
        
        for (dd, cc) in ep :

            x = [ self.approx_middle[DIR.U],
                  self.approx_middle[DIR.U] + dd[DIR.U] * ff ] 
            y = [ self.approx_middle[DIR.V],
                  self.approx_middle[DIR.V] + dd[DIR.V] * ff ] 
            z = [ self.approx_middle[DIR.W],
                  self.approx_middle[DIR.W] + dd[DIR.W] * ff ] 


            ax.plot(x, y, z, zdir='z', color=cc)
        pass


        ## x = [ self.approx_middle[DIR.U],
        ##        self.approx_middle[DIR.U] ]
        ## xp = [ self.approx_middle[DIR.U],
        ##        self.approx_middle[DIR.U] + self. / 4.0 ]
        ## y = [ self.approx_middle[DIR.V],
        ##        self.approx_middle[DIR.V] ]
        ## yp = [ self.approx_middle[DIR.V],
        ##        self.approx_middle[DIR.V] + self. / 4.0 ]
        ## z = [ self.approx_middle[DIR.W],
        ##        self.approx_middle[DIR.W] ]
        ## zp = [ self.approx_middle[DIR.W],
        ##        self.approx_middle[DIR.W] + self.wdir / 4.0  ]

        ## ax.plot(xp, y, z, zdir='z', color='RED')
        ## ax.plot(x, yp, z, zdir='z', color='GREEN')
        ## ax.plot(x, y, zp, zdir='z', color='BLUE')
        
        #ax.annotate('WOW',  

        #plt.annotate('local max', xy=(2, 1), xytext=(3, 1.5),
        #             arrowprops=dict(facecolor='black', shrink=0.05),
        #             ) 


 

    
    pass
    #----------------------------------------------------------------------------------------------------

    def plot_mesh_on_ax(self, ax) :
        
        sys.stdout.flush()

        self.mesh.plot_mesh_on_ax(ax)
    pass


    #----------------------------------------------------------------------------------------------------
    



pass # END OF block CLASS




##########################################################################################################

if __name__ == '__main__' :

    p1 = point(0.0, 0.0, 0.0, {'TEMP': 20, 'VOL_FRAC': 1.0, 'THICK': 1.0 } )
    p2 = point(10.0, 0.0, 0.0, {'TEMP': 40, 'VOL_FRAC': 1.0, 'THICK': 1.0 } )
    p3 = point(10.0, 10.0, 0.0, {'TEMP': 20, 'VOL_FRAC': 1.0, 'THICK': 1.0 } )
    p4 = point(0.0, 10.0, 0.0, {'TEMP': 20, 'VOL_FRAC': 1.0, 'THICK': 1.0 } )


    print 'P1 = ', p1
    print 'P2 = ', p2
    print 'P3 = ', p3
    print 'P4 = ', p4


    v1 = vector(p1, p2)
    v2 = vector(p2, p3)
    v3 = vector(p3, p4)
    v4 = vector(p4, p1)


    #ff = face( 'FACEname', v1, v2, v3, v3 )
    #print 'ff:', str(ff)

    #vv = block('BLOCKname', ff)
    
    #print 'vv:', vv


    print 'v1=', v1
    print 'v2=', v2
    print 'v3=', v3
    print 'v4=', v4

    #print '*' * 80
    #print
    #print 'v1=', v1

    v1f = v1.vector_from_flip()
    #print 'FLIP v1f=', v1f

    v1fs = v1.vector_from_scale(-1.0)
    #print 'SCALE BY-1  v1fs=', v1fs

    #print 'v2=', v2

    nv = v2.cross(v1f).vector_from_length(v2.sp.params['THICK']/2.0)
    #print 'nv = v2 cross v1f =', nv

    nnv = nv.vector_from_flip()
    #print 'nnv = nv flipped', nnv
    
    #print '*' * 80

    fface = face('MID', v1, v2, v3, v4)

    #print 'MID FACE =', fface

    print '=' * 80

    vblock = block('EXV', fface)
    
    vblock.explode_block_from_middle_surface()
    
    print 'VBLOCK = ', vblock

    kl = sorted(vblock.pl.keys())

    for k in kl :
        print k, ' -=>  ', vblock.pl[k]
    pass

#    for k, v in vblock.pl.items() :
#        print k, ' -=>p:  ', v
#    pass
    print '_' * 80

    print 'NUMBER OF EDGES = ', len(vblock.el)
    kl = sorted(vblock.el.keys())
    for k in kl :
        print k, ' -=>  ', vblock.el[k]
    pass
   
    print '=' * 80


    mid_edge_re = re.compile( r'[-]*[1],[ ]*[-]*[1],[ ]*[1],[ ]*[-]*[1],[ ]*[-]*[1],[ ]*[-]*[1]', re.IGNORECASE )


    print 'kl=', kl
    #rkl = filter(mid_edge_re.search, kl)
    fst = filter(lambda f: f[2] == f[5] == -1 , kl)
    #rkl = filter(mid_edge_re.search, kl)


    print 'fst PRESORT=', fst
    fst.sort()
    print 'fst=', fst
 
    for k in fst :
        print k, ' -=>  ', vblock.el[k]
    pass
    


    ## sv = 0.6
    ## lul =  vector(p1, p2)
    ## print 'LUL=', lul
    ## ulop = lul.point_from_scale(sv)
    ## print 'POINT at :', sv, '  -=>', ulop

    ## luh =  vector(p4, p3)
    ## print 'LUH=', luh
    ## uhip = luh.point_from_scale(0.5)
    ## print 'HI HALF WAY POINT:', uhip
    ## tl = vector(ulop, uhip)
    ## print 'TL=', tl
    ## sys.exit(0)
    
    
    lul =  vector(p1, p2)
    lvl =  vector(p1, p4)
    luh =  vector(p4, p3)
    lvh =  vector(p2, p3)

    pp = lul.point_from_length(lul.len*2.0)
    print 'LU.EP=', lul.ep
    print 'EP*2=', pp

    sys.exit(0)



    neu = 3
    nev = 4

    du = 1.0 / float(neu)
    dv = 1.0 / float(nev)

    for u in range(neu+1) :
        uval = du * float(u) 
        ulop = lul.point_from_scale(uval)
        uhip = luh.point_from_scale(uval)
        print 'UUU=', u,  uval
        print 'LOW:', ulop
        print 'HIGH:', uhip

        
        vl = vector(ulop, uhip)
        for v in range(nev+1) :
            vval = dv * float(v)
            pt = vl.scale(vval)
            print 'U=', u, '   V=', v, '  VALS=', pt
        pass
    pass

pass



    
