import sys
import os
import re
import math
import string
import cPickle as pickle
import time
import copy
from datetime import date
from xml.dom import minidom
import numpy as np
from xml.dom import minidom

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm

import constants
from constants import LOC, DIR, DOF, BTAGS, GTAGS, ETAGS, FTAGS, JTAGS, FACE_GRID_INDICES

from collections import namedtuple

from ggen_util import print_timing, enum, here, now, interpolate
from ggen_sup import point, vector, face, block

##########################################################################################################

def write_list_csv(fp, nn, per_line) :
    lnn = len(nn)
    tc = 0
    lc = 0
    for ns in nn :
        tc = tc + 1
        lc = lc + 1
        ns = str(ns)
        ns = ns.strip()
        fp.write(ns)
        if( ( lc == per_line ) or ( tc == lnn ) ) :
            lc = 0
            if( tc != lnn ) :
                fp.write(', ')
            pass
            fp.write('\n')
        else :
            fp.write(', ')
            pass
        pass
    pass
pass

##########################################################################################################

class mat_prop(object) :

    def __init__ (self, mat=None) :
        self.material = mat
        self.temperature = 0.0
        self.props = dict.fromkeys(constants.mat_prop_tags, 0.0)
        self.prop_scale_factor = 1.0      
    pass

    #----------------------------------------------------------------------------------------------------

    def __del__ (self) :
        del self.props
    pass

pass

##########################################################################################################

class material(object) :
    
    def __init__ (self, mname) :
        self.name = mname
        self.lo_temp_defined = None
        self.hi_temp_defined = None
        self.number_of_output_divisions = 5
        # TEMPERATURE DEPENDENT LIST OF MAT PROPS - SORT BASED ON TEMP
        self.mat_props = []
    pass

    #----------------------------------------------------------------------------------------------------

    def __del__ (self) :
        del self.mat_props
    pass

    #----------------------------------------------------------------------------------------------------

    def set_number_of_output_divisions(self, n) :
        self.number_of_output_divisions = max(1, n-1)
    pass

    #----------------------------------------------------------------------------------------------------

    def add_props_at_temp_full(self, temp, density,
                               ym_11, ym_22, ym_33,
                               pr_12, pr_13, pr_23,
                               sm_12, sm_13, sm_23,
                               tec_11, tec_22, tec_33,
                               sy = -1.0, su = -1.0) :
        mp = mat_prop(self)
        mp.temperature = temp
        mp.props['DENSITY'] = density
        mp.props['YM11'] = ym_11
        mp.props['YM22'] = ym_22
        mp.props['YM33'] = ym_33
        mp.props['PR12'] = pr_12
        mp.props['PR13'] = pr_13
        mp.props['PR23'] = pr_23
        mp.props['SM12'] = sm_12
        mp.props['SM13'] = sm_13
        mp.props['SM23'] = sm_23
        mp.props['TEC11'] = tec_11
        mp.props['TEC22'] = tec_22
        mp.props['TEC33'] = tec_33
        mp.props['SY'] = sy
        mp.props['SU'] = su
        self.mat_props.append(mp)
        return(mp)
    pass
        
    #----------------------------------------------------------------------------------------------------

    def add_props_at_temp(self, temp, density, ym, pr, tec, sy = -1.0, su = -1.0) :
        sm = ym / ( 2.0 * ( 1 + pr) ) 
        mp = self.add_props_at_temp_full( temp, density,
                                          ym, ym, ym,
                                          pr, pr, pr,
                                          sm, sm, sm,
                                          tec, tec, tec,
                                          sy, su)
        return(mp)
    pass

    #----------------------------------------------------------------------------------------------------

    def get_mat_props_at_temp_with_volume_fraction_scaling(self, mtemp, volfrac) : 
        # INTERPOLATE/EXTRAPOLATE TO REQUESTED TEMPERATURE
        #     - SCALE YM, DENSITY AND SM WRT volfrac
        #     RETURN  A mat_prop OBJECT
        #  STORE PER ELEMENT / ELSET ????
        
        if( len(self.mat_props) == 0 ) :
            print 'NO PROPERTIES DEFINED FOR MATERIAL: ', self.name
            sys.exit(-1)
        pass

        td_props = copy.deepcopy(self.mat_props[0])
        td_props.temperature = float(mtemp)
        
        #print 'len(self.mat_props)= ', len(self.mat_props)
        
        if( len(self.mat_props) > 1 ) :
            #print 'MORE THAN ONE MAT_PROP DEFINED'
            # EXTRAPOLATE/ INTERPOLATE ON TEMPERATURE FIRST
            ts = []
            for tdp in self.mat_props :
                ts.append(tdp.temperature)
            pass
            rts = list(reversed(ts))
            for pv in constants.mat_prop_tags :
                vs = []
                for tdp in self.mat_props :
                    vs.append(tdp.props[pv])
                pass
                #print '\nPV=', pv
                #print 'mtemp =', mtemp
                #print 'vs=', vs
                #print 'ts=', ts
                
                loval = interpolate(vs[0], ts[0], vs[1], ts[1], float(mtemp))

                rvs = list(reversed(vs))
                hival = interpolate(rvs[0], rts[0], rvs[1], rts[1], float(mtemp))
                
                #print 'rvs=', rvs
                #print 'rts=', rts

                #print 'loval = ;', loval
                #print 'hival = ;', hival

                tv = np.interp(mtemp, ts, vs, loval, hival)
                #print 'tv =', tv
                td_props.props[pv] = tv
            pass
        pass                   


        # SCALE DENSITY AND MODULI WRT VOLUME FRACTION
        td_props.props['DENSITY'] *= float(volfrac)
        td_props.props['YM11'] *= float(volfrac)
        td_props.props['YM22'] *= float(volfrac)
        td_props.props['YM33'] *= float(volfrac)
        td_props.props['SM12'] *= float(volfrac)
        td_props.props['SM13'] *= float(volfrac)
        td_props.props['SM23'] *= float(volfrac)

        td_props.prop_scale_factor = float(volfrac)


        return(td_props)
   
        pass
    pass


    #----------------------------------------------------------------------------------------------------

pass

##########################################################################################################
# GRID IS DERIVED FROM POINT

class grid(point) :

    def __init__ (self, gid = 0, x = 0.0, y = 0.0 , z = 0.0, rest = {}) :
        super(grid, self).__init__(x, y, z, rest)
        self.name = 'G' + str(gid)
        self.id = gid
        self.bc = set()
        self.norm_uvw = [0.0, 0.0, 0.0]
        self.uvw = [0, 0, 0]
        
        self.mesh = None
        self.el = set()
        self.results = {}
        self.max_results = {}
        self.min_results = {}
        self.avg_results = {}
    pass

    #----------------------------------------------------------------------------------------------------

    def __del__ (self) :
        super(grid, self).__del__()
        del self.bc
        del self.el
        del self.results
        del self.max_results
        del self.min_results
        del self.avg_results
    pass

    #----------------------------------------------------------------------------------------------------

    def calc_avg_parameter_data_from_gl(self, gl) :
        #print '>>> TOP >>> ', here(), '  NAME = ',self.name
        
        super(grid, self).__init__(0.0, 0.0, 0.0)

        self.calc_avg_parameter_data_from_point_list(gl)

        ## ng = len(gl)
        ## if( ng == 0 ) : return

        
        ## #print 'GRID:', self

        ## #print 'GL = ', gl
        ## # AVERAGE X, Y, Z
        ## for g in gl :
        ##     #print 'TYPE self = ', type(self), '   TYPE g =', type(g)
        ##     for dd in range(len(self.v)) :
        ##         self.v[dd] += g.v[dd]
        ##     pass
        ## pass
        ## for dd in range(len(self.v)) :
        ##      self.v[dd] /= float(ng) 
        ## pass

        ## # AVERAGE ANCILLARY DATA - IF POSSIBLE
        ## pd = {}

        ## # SUM
        ## for g in gl :
        ##     #print 'LOOKING AT GRID PARAMS FOR GRID;', g
        ##     for k, v in g.params.items() :
        ##         #print '  LOOKING AT PARAM :', k
        ##         if( k not in pd ) :
        ##             pd[k] = [v, 1]  # [ sum, count ] pairs
        ##             #print 'FIRST TIME FOR PARAM  VAL=', v
        ##         else :
        ##             if( not isinstance(v, str) ) : # CAN'T SUM STRING PARAMS
        ##                 pd[k][0] += v
        ##                 pd[k][1] += 1
        ##                 #print pd[k][1], ' TIME FOR PARAM  VAL=', v, '   RUNNING TOT = ', pd[k][0] 
        ##             pass
        ##         pass
        ##     pass
        ## pass

        ## # NOW AVERAGE
        ## for k, v  in pd.items() :
        ##     if( not isinstance(v[0], str) ) :
        ##         self.set_param(k, v[0] / float(v[1]))
        ##     else :
        ##        self.set_param(k, v[0])
        ##     pass
        ##     #print 'PARAM ', k, '  AVG VAL = ', self.get_param(k, -1.0)
        ## pass

    pass
        


    #----------------------------------------------------------------------------------------------------

    def __repr__(self) :
        td = GTAGS._asdict()
        ts = ''
        for ( k, v ) in td.items() :
            if( v & self.tags ) :
                if( len(ts) > 0 ) : ts = ts + ','
                ts = ts + k   
            pass
        pass
        tsl = ts.split(',')
        tsl.sort()
        ts = ', '.join(tsl)

        bc_tags = []    
        bcl = sorted(list(self.bc))
        for dof in range(len(constants.dof_tags)) :
            bc_tags.append(0)
        pass
        for dof in bcl :
            bc_tags[dof-1] = 1
        pass

        bc_stags = ''
        for dof in range(len(constants.dof_tags)):
            bc_stags = bc_stags + str(bc_tags[dof])
        pass
        

        mn = ' '
        if( self.mesh is not None ) :
            mn = self.mesh.name
        pass
        
        
#        f = 'G: {0} ->  {1}, {2}, {3}    ( {4}, {5} ) {6:0' + str(len(td)) + 'b} ({7})'
#        s = f.format(self.id, self.v[LOC.X], self.v[LOC.Y], self.v[LOC.Z],
#                     self.uvw[DIR.U], self.uvw[DIR.V], self.tags, ts)
        f = 'G: {0} ->  {1}, {2}, {3}    ( {4}, {5}, {6} ) BC:{7} ({8}) {9}'
        s = f.format(self.id, self.v[LOC.X], self.v[LOC.Y], self.v[LOC.Z],
                     self.uvw[DIR.U], self.uvw[DIR.V], self.uvw[DIR.W], bc_stags, ts, mn)
        return(s)
    pass

    #----------------------------------------------------------------------------------------------------

    def save_abq(self, fp) :
        if( self.id <= 0 ) : return
        if( self.tags & GTAGS.REMOVED ) : return
        f = '{0}, {1}, {2}, {3}'
        s = f.format(self.id, self.v[LOC.X], self.v[LOC.Y], self.v[LOC.Z])
        fp.write(s + '\n')
    pass
        
    #----------------------------------------------------------------------------------------------------


pass

##########################################################################################################

class element(object) :
    def __init__ (self, eid) :
        self.name = 'E' + str(eid)
        self.single_element_set_name = self.name
        self.id = eid
        self.gl = []
        self.tags = 0
        self.results = {}
        self.mat_name = ""
        self.uvw = [0, 0, 0]
        self.aspect_ratio = -1.0
        self.avg_grid = grid( (-1) * self.id ) # TAG THE avg_grid WITH NEG ELEMENT ID
        self.mesh = None
    pass

    #----------------------------------------------------------------------------------------------------

    def __del__ (self) :
        del self.gl
        del self.results
        del self.avg_grid
    pass

    #----------------------------------------------------------------------------------------------------

    def __repr__(self) :
        b2b_brace_re = re.compile( r'[}][ ]*[{]', re.IGNORECASE )

        td = FTAGS._asdict()
        ts = ''
        for ( k, v ) in td.items() :
            if( v & self.tags ) :
                if( len(ts) > 0 ) : ts = ts + ','
                ts = ts + k   
            pass
        pass
        tsl = ts.split(',')
        tsl.sort()
        ts = ', '.join(tsl)
        
        ss = '{:6}' * len(self.gl)
        ff = 'E: |{:6}| ' + '({})' + ss 
        # PUT COMMAS BETWEEN THE FIELDS
        f =  b2b_brace_re.sub(r'}, {', ff)
        gids = []
        for g in self.gl :
            #print 'g TYPE = ', type(g)
            gids.append(g.id)
        pass
        fs = f.format(self.id, ts, *gids)
        #fs = f.format(self.id, *gids) + '\n'
        #print s
        width = 80
        fields_per_line = 3
        field_width = int( width / fields_per_line )
        #print 'self.results keys= ', self.results.keys()
        for ip, res in self.results.items() :
            cnt = 0
            tcnt = 0
            fs += 'INTEGRATION POINT: ' + str(ip) + '\n'
            #print 'res =', res
            rl = len(res)
            kl = sorted(res.keys())
            for k in kl :
                v = res[k]
                #print 'KEY:', k, '   VAL:', v
                s = k + ' = ' + str(v)
                s = s.rjust(field_width)
                fs += fs + '\n' + s 
                cnt = cnt + 1
                tcnt = tcnt + 1

                if( ( cnt == fields_per_line ) or ( tcnt == rl ) ) :
                    #print
                    fs += fs + '\n' 
                    cnt = 0
                pass
            pass
        pass
        return(fs)
    pass

    #----------------------------------------------------------------------------------------------------

    def save_abq(self, fp) :
        if( self.id <= 0 ) : return
        
        fp.write(str(self.id) + ', ')
        tc = 0
        lc = 1
        for g in self.gl :
            tc = tc + 1
            lc = lc + 1
            fp.write(str(g.id))
            if( ( lc == 13 ) or ( tc == len(self.gl) ) ) :
                lc = 0
                if( tc != len(self.gl) ) :
                    fp.write(', ')
                pass
                fp.write('\n')
            else :
                if( tc != len(self.gl) ) :
                    fp.write(', ')
                    pass
                pass
            pass
        pass
    pass

    #----------------------------------------------------------------------------------------------------

pass


##########################################################################################################

class elset(object) :

    def __init__ (self, iname, eindex, mm) :
        self.name = iname
        self.index = eindex
        self.mesh = mm
        self.el = []
        self.avg_grid = grid()
    pass

    #----------------------------------------------------------------------------------------------------

    def __del__ (self) :
        del self.el
        del self.avg_grid
    pass

    #----------------------------------------------------------------------------------------------------

    def calc_elset_parameter_averages(self) :
        gl = []
        for e in self.el :
            gl.append(e.avg_grid)
        pass
        self.avg_grid.calc_avg_parameter_data_from_gl(gl)
    pass
    
    #----------------------------------------------------------------------------------------------------

pass
    
##########################################################################################################


class mesh(object) :

    def __init__ (self, mname) :
        self.name = mname
        self.eset_name = 'E'+self.name
        self.nset_name = 'N'+self.name
        
        self.model = None
        self.block = None
        
        self.el = {}
        self.gl = {}

        self.u_elsets = []
        self.v_elsets = []
        self.w_elsets = []
        
        self.material = None
        self.temp_lo = None
        self.temp_hi = None
    pass

    #----------------------------------------------------------------------------------------------------

    def __del__ (self) :
        del self.el
        del self.gl
        del self.u_elsets
        del self.v_elsets
        del self.material
    pass

    #----------------------------------------------------------------------------------------------------
    # http://stackoverflow.com/questions/1259963/python-assert-that-variable-is-instance-method/1260881#1260881

    ## def is_instance_method(self, obj):
    ##     """Checks if an object is a bound method on an instance."""
    ##     if( not isinstance(obj, types.MethodType) ) :
    ##         return(False) # Not a method
    ##     pass
    ##     if( obj.im_self is None ) :
    ##         return(False) # Method is not bound
    ##     pass
    ##     if( ( issubclass(obj.im_class, type) )
    ##         or ( obj.im_class is types.ClassType ) ) :
    ##         return(False) # Method is a classmethod
    ##     pass
    ##     return(True)
    ## pass

    ## #----------------------------------------------------------------------------------------------------

    ## def __getstate__(self):
    ##     return( dict((k, v) for k, v in self.__dict__.iteritems()
    ##                    if not self.is_instance_method(getattr(self, k))) )
    ## pass

    #----------------------------------------------------------------------------------------------------

    def set_material(self, mat) :
        self.material = mat
    pass

                 
    #----------------------------------------------------------------------------------------------------

    def plot_mesh_on_ax(self, ax) :
        #TBD
        pass
    pass

     

    #----------------------------------------------------------------------------------------------------

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


    def generate(self) :
        n = -1
        p = 1
        
        b = self.block
        b.ng = b.ne * 2 + 1
        #print 'NG =', b.ng
        #print 'NE =', b.ne



        center_line_indices = [[], [], []]

        vl = [None] * 4 # WORKING VECTOR LIST
        pl = [None] * 7 # WORKING POINT LIST

        duvw = np.zeros(3)
        
        for i in range(3) :
            cuvw = b.ng[i] / 2 
            center_line_indices[i].append(cuvw)
            if( ( b.ng[i] % 2 ) == 0 ) : # EVEN NUMBER OF GRIDS THEN TAG 2 ROWS AS CENTER GRIDS
                center_line_indices[i].append( cuvw + 1 )
            pass
            duvw[i] = 1.0 / float(b.ng[i] - 1.0)        
        pass

        uvw = np.zeros(3)
        iuvw = np.zeros(3, dtype=np.int32)
        for iu in range(b.ng[DIR.U]) :
            
            iuvw[DIR.U] = iu
            uvw[DIR.U] = duvw[DIR.U] * float(iu)
            # AVOID ROUND OFF
            if( iu == 0 ) : uvw[DIR.U] = 0.0
            if( iu == ( b.ng[DIR.U] - 1 ) ) : uvw[DIR.U] = 1.0
            
            # SCALE THE U DIRECTION EDGE VECTORS
            pl[0] = b.el[n, n, n, p, n, n].point_from_scale(uvw[DIR.U])
            pl[1] = b.el[n, p, n, p, p, n].point_from_scale(uvw[DIR.U])
            pl[2] = b.el[n, p, p, p, p, p].point_from_scale(uvw[DIR.U])
            pl[3] = b.el[n, n, p, p, n, p].point_from_scale(uvw[DIR.U])

            # VECTORS ON THE +/- VW CUBE FACES POINTING IN THE POSITIVE V DIRECTION
            vl[0] = vector(pl[0], pl[1])
            vl[1] = vector(pl[3], pl[2])
            
            for iv in range(b.ng[DIR.V]) :
                iuvw[DIR.V] = iv
                uvw[DIR.V] = duvw[DIR.V] * float(iv)
                # AVOID ROUND OFF
                if( iv == 0 ) : uvw[DIR.V] = 0.0
                if( iv == ( b.ng[DIR.V] - 1 ) ) : uvw[DIR.V] = 1.0
                
                # SCALE THE V DIRECTION VECTORS
                pl[4] = vl[0].point_from_scale(uvw[DIR.V])
                pl[5] = vl[1].point_from_scale(uvw[DIR.V])

                # VECTOR ACROSS THE VW CUT PLANE IN THE POSITIVE W DIRECTION
                vl[3] = vector(pl[4], pl[5])
                
                for iw in range(b.ng[DIR.W]) :


                    # BEGIN FIGURE TAGS 
                    gtags = 0

                    if( iu == 0 ) : gtags |= GTAGS.MIN_U
                    if( iu == (b.ng[DIR.U]-1) ) : gtags |= GTAGS.MAX_U
                    if( iu in center_line_indices[DIR.U] ) : gtags |= GTAGS.CENTER_U

                    if( iv == 0 ) : gtags |= GTAGS.MIN_V
                    if( iv == (b.ng[DIR.V]-1) ) : gtags |= GTAGS.MAX_V
                    if( iv in center_line_indices[DIR.V] ) : gtags |= GTAGS.CENTER_V

                    if( iw == 0 ) : gtags |= GTAGS.MIN_W
                    if( iw == (b.ng[DIR.W]-1) ) : gtags |= GTAGS.MAX_W
                    if( iw in center_line_indices[DIR.W] ) : gtags |= GTAGS.CENTER_W

                    # TAG MID GRIDS
                    if( ( iu % 2 == 0 ) and ( iv % 2 == 0 ) and ( iw % 2 == 0 ) ) :
                        gtags |= GTAGS.ELEMENT_CORNER
                    pass
                    if( ( iu % 2 == 1 ) and ( iv % 2 == 1 ) and ( iw % 2 == 1 ) ) :
                        gtags |= GTAGS.ELEMENT_MID_VOL
                        gtags |= GTAGS.ANCILLARY
                        #gtags |= GTAGS.PENDING_DELETE
                    pass
                        
                    if( ( iu % 2 == 0 ) and ( iv % 2 == 0 ) and ( iw % 2 == 1 ) ) :
                        gtags |= GTAGS.ELEMENT_MID_EDGE
                        gtags |= GTAGS.ELEMENT_MID_EDGE_W
                    pass
                    if( ( iu % 2 == 0 ) and ( iv % 2 == 1 ) and ( iw % 2 == 1 ) ) :
                        gtags |= GTAGS.ELEMENT_MID_FACE
                        gtags |= GTAGS.ELEMENT_MID_FACE_VW
                        gtags |= GTAGS.PENDING_DELETE
                    pass
                        
                    if( ( iu % 2 == 0 ) and ( iv % 2 == 1 ) and ( iw % 2 == 0 ) ) :
                        gtags |= GTAGS.ELEMENT_MID_EDGE
                        gtags |= GTAGS.ELEMENT_MID_EDGE_V
                    pass
                    if( ( iu % 2 == 1 ) and ( iv % 2 == 1 ) and ( iw % 2 == 0 ) ) :
                        gtags |= GTAGS.ELEMENT_MID_FACE
                        gtags |= GTAGS.ELEMENT_MID_FACE_UV
                        gtags |= GTAGS.PENDING_DELETE
                    pass
                        
                    if( ( iu % 2 == 1 ) and ( iv % 2 == 0 ) and ( iw % 2 == 0 ) ):
                        gtags |= GTAGS.ELEMENT_MID_EDGE
                        gtags |= GTAGS.ELEMENT_MID_EDGE_U
                    pass
                    if( ( iu % 2 == 1 ) and ( iv % 2 == 0 ) and ( iw % 2 == 1 ) ) :
                        gtags |= GTAGS.ELEMENT_MID_FACE
                        gtags |= GTAGS.ELEMENT_MID_FACE_UW
                        gtags |= GTAGS.PENDING_DELETE
                    pass


                    # DON'T DEFINE THE GRIDS WE DON'T NEED
                    if( gtags & GTAGS.PENDING_DELETE ) :
                        continue
                    pass
                
                    # END FIGURE TAGS


                    iuvw[DIR.W] = iw
                    uvw[DIR.W] = duvw[DIR.W] * float(iw)
                    # AVOID ROUND OFF
                    if( iw == 0 ) : uvw[DIR.W] = 0.0
                    if( iw == ( b.ng[DIR.W] - 1 ) ) : uvw[DIR.W] = 1.0
                    
                    # SCALE THE W DIRECTION VECTOR
                    pt = vl[3].point_from_scale(uvw[DIR.W])

                    self.model.gid_max += 1
                    gid = self.model.gid_max
                

                    g = grid(gid, pt.v[LOC.X], pt.v[LOC.Y], pt.v[LOC.Z], pt.params)
                    for d in [ DIR.U, DIR.V, DIR.W ] :
                        g.norm_uvw[d] = uvw[d]
                    pass
                    #print 'NORM_UVW = ', g.norm_uvw
                    g.uvw[DIR.U] = iu
                    g.uvw[DIR.V] = iv
                    g.uvw[DIR.W] = iw
                    g.tags = gtags
                    g.mesh = self
                    self.gl[iu, iv, iw] = g
                pass # W LOOP
            pass # V LOOP
        pass # U LOOP


        #print 'GOING TO DEFINE SOME ELSETS FOR EACH INCREMENT IN EACH DIR'

        self.u_elsets = []
        for iu in range(b.ne[DIR.U]) :
            s = self.name + 'EU' + str(iu)
            els = elset(s, iu, self)
            self.u_elsets.append(els)
        pass
    
        self.v_elsets = []
        for iv in range(b.ne[DIR.V]) :
            s = self.name + 'EV' + str(iv)
            els = elset(s, iv, self)
            self.v_elsets.append(els)
        pass

        self.w_elsets = []
        for iw in range(b.ne[DIR.W]) :
            s = self.name + 'EW' + str(iw)
            els = elset(s, iw, self)
            self.w_elsets.append(els)
        pass

        #print 'GOING TO DEFINE THE ELEMENTS'
        
        for iu in range(b.ne[DIR.U]) :
            for iv in range(b.ne[DIR.V]) :
                for iw in range(b.ne[DIR.W]) :

                    etags = 0

                    # FROM CALCULIX MANUAL...
                    #  face 1: 1-2-3-4 -> FTAGS.MIN_W - 1
                    #  face 2: 5-8-7-6 -> FTAGS.MAX_W - 2
                    #  face 3: 1-5-6-2 -> FTAGS.MIN_V - 4
                    #  face 4: 2-6-7-3 -> FTAGS.MAX_U - 8
                    #  face 5: 3-7-8-4 -> FTAGS.MAX_V - 16
                    #  face 6: 4-8-5-1 -> FTAGS.MIN_U - 32

                    #face_num = (math.log(etag) / math.log(2.0)) + 1

                    # TAG WHICH ELEMENTS ARE ON THE FACES OF THE BLOCK
                    if( iu == 0 ) : etags |= FTAGS.MIN_U
                    if( iu == (b.ne[DIR.U]-1) ) : etags |= FTAGS.MAX_U

                    if( iv == 0 ) : etags |= FTAGS.MIN_V
                    if( iv == (b.ne[DIR.V]-1) ) : etags |= FTAGS.MAX_V

                    if( iw == 0 ) : etags |= FTAGS.MIN_W
                    if( iw == (b.ne[DIR.W]-1) ) : etags |= FTAGS.MAX_W


                    self.model.eid_max += 1
                    eid = self.model.eid_max
                    e = element(eid)
                    e.tags = etags
                    
                    # GRID INDICES
                    iug = iu * 2
                    ivg = iv * 2
                    iwg = iw * 2

                    e.gl.append(self.gl[iug, ivg, iwg])      
                    e.gl.append(self.gl[iug+2, ivg, iwg])    
                    e.gl.append(self.gl[iug+2, ivg+2, iwg])
                    e.gl.append(self.gl[iug, ivg+2, iwg])

                    e.gl.append(self.gl[iug, ivg, iwg+2])
                    e.gl.append(self.gl[iug+2, ivg, iwg+2])
                    e.gl.append(self.gl[iug+2, ivg+2, iwg+2])
                    e.gl.append(self.gl[iug, ivg+2, iwg+2])

                    e.gl.append(self.gl[iug+1, ivg, iwg])
                    e.gl.append(self.gl[iug+2, ivg+1, iwg])
                    e.gl.append(self.gl[iug+1, ivg+2, iwg])
                    e.gl.append(self.gl[iug, ivg+1, iwg])

                    e.gl.append(self.gl[iug+1, ivg, iwg+2])
                    e.gl.append(self.gl[iug+2, ivg+1, iwg+2])
                    e.gl.append(self.gl[iug+1, ivg+2, iwg+2])
                    e.gl.append(self.gl[iug, ivg+1, iwg+2])

                    e.gl.append(self.gl[iug, ivg, iwg+1])      
                    e.gl.append(self.gl[iug+2, ivg, iwg+1])    
                    e.gl.append(self.gl[iug+2, ivg+2, iwg+1])
                    e.gl.append(self.gl[iug, ivg+2, iwg+1])

                    e.uvw[DIR.U] = iu
                    e.uvw[DIR.V] = iv
                    e.uvw[DIR.W] = iw
                    e.mesh = self


                    # UPDATE THE BACKPOINTER LIST OF ELEMENTS THAT EACH GRID BELONGS TO
                    for g in e.gl :
                        g.el.add(e)
                    pass


                    self.el[iu, iv, iw] = e

                    self.u_elsets[iu].el.append(e)
                    self.v_elsets[iv].el.append(e)
                    self.w_elsets[iw].el.append(e)

                    


                pass # W LOOP
            pass # V LOOP
        pass # U LOOP

        #print 'DONE DEFINING ELEMENTS'
    
        #self.print_grids()
        #self.print_elements()
    
    pass


    #----------------------------------------------------------------------------------------------------

    def print_grids(self) :
        
        kl = sorted(self.gl.keys())

        cnt = 0
        
        print '_'*80 + '\n'
        print 'GRIDS: MESH: |' + self.name + '| # OF SIZE OF MESH CLASS GRID DICT= ', len(self.gl)
        for k in kl :
            if( self.gl[k] is None ) : continue
            print k, '-', cnt, ' -=>  ', self.gl[k]
            cnt += 1
        pass
        print '# OF POPULATED ELEMENTS IN MESH CLASS GRID DICT = ', cnt, ' (', float(cnt)/float(len(self.gl)) * 100.0, '% )'
        print '_'*80 + '\n'
        sys.stdout.flush()
    pass

    #----------------------------------------------------------------------------------------------------

    def print_elements(self) :
        
        kl = sorted(self.el.keys())

        cnt = 0
        print '_'*80 + '\n'
        print 'ELEMS: MESH: |' + self.name + '| # OF SIZE OF MESH CLASS ELEMENT DICT= ', len(self.el)
        for k in kl :
            if( self.el[k] is None ) : continue
            print k, '-', cnt, ' -=>  ', self.el[k]
            cnt += 1
        pass
        print '# OF POPULATED ELEMENTS IN MESH CLASS ELEMENT DICT = ', cnt, ' (', float(cnt)/float(len(self.el)) * 100.0, '% )'
        print '_'*80 + '\n'
        sys.stdout.flush()
    pass

    #----------------------------------------------------------------------------------------------------

    def get_grid_list_from_tags(self, tag, igl = None) :
        # If a grid list is passed in only look in that list - otherwise consider all the grids
        rgl = self.gl.values()
        if( igl is not None ) :
            rgl = igl
        pass
        rgl = filter( ( lambda gg: gg is not None ), rgl)
        rgl = filter( ( lambda gg: gg.tags & tag ), rgl )
        return(rgl)
    pass

    #----------------------------------------------------------------------------------------------------

    def get_grid_list_from_not_tags(self, tag, igl = None) :
        # If a grid list is passed in only look in that list - otherwise consider all the grids
        rgl = self.gl.values()
        if( igl is not None ) :
            rgl = igl
        pass
        rgl = filter( ( lambda gg: gg is not None ), rgl)
        rgl = filter( ( lambda gg: not ( gg.tags & tag ) ), rgl )
        return(rgl)
    pass

    #----------------------------------------------------------------------------------------------------

    def get_grid_from_indices(self, iu, iv, iw = 0, igl = None) :
        #print 'LOOKING FOR GRID AT INDEX =', iu, iv, iw

        rgl = self.gl;
        if( igl is not None ) :
            rgl = igl
        pass

        
        rgl = filter( ( lambda gg: (gg.uvw[DIR.U] == iu) and ( gg.uvw[DIR.V] == iv ) and ( gg.uvw[DIR.W] == iw ) ), rgl )

        rg = None
        if( len(rgl) > 0 ) :
            rg = rgl[0]
        pass
        #print len(rgl), rgl, rg
        #if( rg is None ) :
        #    print 'ERROR'
        #    sys.exit(-1)
        #pass
        return(rg)
    pass

    #----------------------------------------------------------------------------------------------------

    def get_grid_list_from_indices(self, iu = None, iv = None, iw = None, igl = None) :
        #print '>>> TOP >>> ', here(), '  NAME = ',self.name
        #  REF: PYTHON library reference
        # Shallow copies of dictionaries can be made using dict.copy(),
        #   and of lists by assigning a slice of the entire list, for
        #   example, copied_list = original_list[:]        

        rgl = self.gl;
        if( igl is not None ) :
            rgl = igl
        pass

    
        #print 'LOOKING FOR GRID(S) AT INDEX =', iu, iv, iw
        #print 'BEFORE RGL = '
        #for g in rgl :
        #   print g.uvw[DIR.U], g.uvw[DIR.V], g.uvw[DIR.W], '->', g
        #pass

    
        if( iu is not None ) :
            rgl = filter( ( lambda gg: gg.uvw[DIR.U] == iu ), rgl )
        pass
        if( iv is not None ) :
            rgl = filter( ( lambda gg: gg.uvw[DIR.V] == iv ), rgl )
        pass
        if( iw is not None ) :
            rgl = filter( ( lambda gg: gg.uvw[DIR.W] == iw ), rgl )
        pass

        ## print 'AFTER RGL = '
        ## for g in rgl :
        ##    print g.uvw[DIR.U], g.uvw[DIR.V], g.uvw[DIR.W], '->', g
        ## pass

    
        return(rgl)
    pass

    #----------------------------------------------------------------------------------------------------

    def get_element_list_from_indices(self, iu = None, iv = None, iw = None) :
        #print 'LOOKING FOR ELEMENT AT INDEX =', iu, iv, iw
        #  REF: PYTHON library reference
        # Shallow copies of dictionaries can be made using dict.copy(),
        #   and of lists by assigning a slice of the entire list, for
        #   example, copied_list = original_list[:]        
        rel = self.el[:];
        if( iu is not None ) :
            rel = filter( ( lambda ee: ee.uvw[DIR.U] == iu ), rel )
        pass
        if( iv is not None ) :
            rel = filter( ( lambda ee: ee.uvw[DIR.V] == iv ), rel )
        pass
        if( iw is not None ) :
            rel = filter( ( lambda ee: ee.uvw[DIR.W] == iw ), rel )
        pass
        return(rel)
    pass

    #----------------------------------------------------------------------------------------------------

    def get_elements_using_grid_id(self, gid) :
        gels = set() # USE set TO MAKE SURE WE BUILD A UNIQUE LIST
        for ee in self.el :
            for gg in ee.gl :
                if( gg.id == gid ) :
                    gels.add(ee)
                    break
                pass
            pass
        pass
        gel = list(gels)
        return(gel)
    pass

    #----------------------------------------------------------------------------------------------------

    def get_elements_using_grid(self, g) :
        gels = set() # USE set TO MAKE SURE WE BUILD A UNIQUE LIST
        for ee in self.el :
            for gg in ee.gl :
                if( gg == g ) :
                    gels.add(ee)
                    break
                pass
            pass
        pass
        gel = list(gels)
        return(gel)
    pass

    #----------------------------------------------------------------------------------------------------

    def calc_element_parameter_averages(self) :
        #print '>>> TOP >>> ', here()
        
        self.temp_lo = self.model.initial_temp
        self.temp_hi = self.model.initial_temp
        for e in self.el.values() :
            #print 'LOOKING AT ELEMENT ', e.id, '   LEN ELEMENT GL = ', len(e.gl)
            e.avg_grid.calc_avg_parameter_data_from_gl(e.gl)
            for g in e.gl :
                temp = g.get_param('TEMP', self.model.initial_temp)
                #if( temp is None ) : continue
                #print '  LOOKING AT GRID ', g.id, ' T=', g.temp 
                self.temp_lo = min(self.temp_lo, temp)
                self.temp_hi = max(self.temp_hi, temp)
                #print '   ', self.temp_lo, '  ', self.temp_hi
            pass
        pass

        # AVERAGE GRID DATA FOR ELSET (EACH ELEMENT ROW)
        for es in ( self.u_elsets + self.v_elsets + self.w_elsets ):
            es.calc_elset_parameter_averages()
        pass
    pass


    #----------------------------------------------------------------------------------------------------

    def save_abq_grids(self, fp) :
        tgl = self.gl.values()
        tgl = filter( ( lambda g: g is not None ), tgl)
        tgl.sort(key=lambda gg : gg.id)
        for gg in tgl :
            if( gg.id <= 0 ) : continue
            gg.save_abq(fp)
        pass
    pass

    #----------------------------------------------------------------------------------------------------
    # TO BE COMPLETELY GENERAL WE WILL HAVE TO MAKE A ELSET FOR EVERY ELEMENT
    
    def save_abq_elements(self, fp) :
        #elsets_list = []
        tel = self.el.values()
        tel.sort(key=lambda ee : ee.id)
        
        for ee in tel :
            fp.write('*ELEMENT, TYPE=C3D20R, ELSET=' + ee.single_element_set_name + '\n')
            ee.save_abq(fp)
        pass
    pass

    #----------------------------------------------------------------------------------------------------

pass


##########################################################################################################

class join_face(object) :

    #----------------------------------------------------------------------------------------------------

    def __init__(self, ref_block, ref_block_tag, parent_join) :
        self.parent = parent_join
        self.tags = 0
        self.block = ref_block
        self.face_grid_tag = ref_block_tag
        self.hinge_dir = None

        self.face_el = set()

        
        # GRID LIST OF GRIDS ON HINGE LINE - FOUND NO MATTER WHAT KIND OF JOIN IT IS
        self.gl = []
        self.ortho_gl = []
        
        # GRID LIST OF GRIDS ON JOIN FACE - FOUND NO MATTER WHAT KIND OF JOIN IT IS
        self.face_gl = []

        self.face_dirs = [-1, -1]

        self.get_face_directions()
        
        self.name = self.parent.name + '_' + self.block.name
    pass

    #----------------------------------------------------------------------------------------------------
    # THIS IS CALLED BY __init__
    def get_face_directions(self) :
        self.face_dirs = [-1, -1]
        if( ( self.face_grid_tag == GTAGS.MAX_U )
            or ( self.face_grid_tag == GTAGS.MIN_U ) ) :
            self.face_dirs = [ DIR.V, DIR.W ]
        pass
        if( (  self.face_grid_tag == GTAGS.MAX_V )
            or (  self.face_grid_tag == GTAGS.MIN_V ) ) :
             self.face_dirs = [ DIR.U, DIR.W ]
        pass
        if( (  self.face_grid_tag == GTAGS.MAX_W )
            or (  self.face_grid_tag == GTAGS.MIN_W ) ) :
             self.face_dirs = [ DIR.U, DIR.V ]
        pass
    pass

    #----------------------------------------------------------------------------------------------------
    # WE REALLY DON'T NEED THIS
    def get_face_lines(self) :
        if( self.face_dirs[1] == DIR.U ) :
            main_grid_line_tag = GTAGS.MIN_U
        pass
        if( self.face_dirs[1] == DIR.V ) :
            main_grid_line_tag = GTAGS.MIN_V
        pass
        if( self.face_dirs[1] == DIR.W ) :
            main_grid_line_tag = GTAGS.MIN_W
        pass
        if( self.face_dirs[0] == DIR.U ) :
            ortho_grid_line_tag = GTAGS.MIN_U
        pass
        if( self.face_dirs[0] == DIR.V ) :
            ortho_grid_line_tag = GTAGS.MIN_V
        pass
        if( self.face_dirs[0] == DIR.W ) :
            ortho_grid_line_tag = GTAGS.MIN_W
        pass

        ma = self.block.mesh
        self.gl = ma.get_grid_list_from_tags(main_grid_line_tag,
                                             self.face_gl)
            
        self.ortho_gl =  ma.get_grid_list_from_tags(ortho_grid_line_tag,
                                                    self.face_gl)
    pass

    #----------------------------------------------------------------------------------------------------

    def __del__ (self) :
        del self.gl
        del self.face_gl
    pass

    #----------------------------------------------------------------------------------------------------
        
pass


##########################################################################################################

class join(object) :

    #----------------------------------------------------------------------------------------------------

    def __init__(self, jname, ba, ba_tag, bb, bb_tag, jtype, parent_model) :
        self.name = jname
        self.jface_a = join_face(ba, ba_tag, self)
        self.jface_b = join_face(bb, bb_tag, self)
        self.jfaces = [ self.jface_a, self.jface_b ]
        self.type = jtype

        self.slave_tag = ''

        self.model = parent_model

        # THESE ARE SET EQUAL TO a_gl OR b_gl DEPENDING ON THIER LENGTH (slave_gl IS THE LONGER)
        self.slave = None 
        self.master = None
        
        self.edge_vector = np.zeros(3)
        
    pass

    #----------------------------------------------------------------------------------------------------

    def __del__ (self) :
        del self.jface_a
        del self.jface_b
    pass

    #----------------------------------------------------------------------------------------------------

    def set_a_hinge_direction(self, dir_tag) :
        self.jface_a.hinge_dir = dir_tag
    pass

    #----------------------------------------------------------------------------------------------------

    def set_b_hinge_direction(self, dir_tag) :
        self.jface_b.hinge_dir = dir_tag
    pass



    #----------------------------------------------------------------------------------------------------
    # STATIC
    # GIVEN THE BLOCK FACE TAG AND THE HINGE DIRECTION RETURN THE TAG THAT WILL RETRIEVE THE HINGE LINE NODES
    @staticmethod
    def get_face_hinge_line_tags(block_face_grid_tag, hinge_dir_tag) :

        hinge_grid_tag = -1
        ortho_grid_tag = -1
        ortho_grid_select_tag = -1

        # TBD - SHOULD SEE IF WE CAN "AND" THE ORTHO GRID TAGS  - SHOULD SPEED THINGS US 
        
        if( block_face_grid_tag == GTAGS.MIN_U  ) :
            if( hinge_dir_tag == DIR.V ) :
                hinge_grid_tag = GTAGS.CENTER_W
                ortho_grid_tag = GTAGS.MIN_V
                ortho_grid_select_tag = GTAGS.MIN_W
            pass
            if( hinge_dir_tag == DIR.W ) :
                hinge_grid_tag = GTAGS.CENTER_V
                ortho_grid_tag = GTAGS.MIN_W
                ortho_grid_select_tag = GTAGS.MAX_V
            pass
        pass
    
        if( block_face_grid_tag == GTAGS.MAX_U ) :
            if( hinge_dir_tag == DIR.V ) :
                hinge_grid_tag = GTAGS.CENTER_W
                ortho_grid_tag = GTAGS.MIN_V
                ortho_grid_select_tag = GTAGS.MAX_W
            pass
            if( hinge_dir_tag == DIR.W ) :
                hinge_grid_tag = GTAGS.CENTER_V
                ortho_grid_tag = GTAGS.MIN_W
                ortho_grid_select_tag = GTAGS.MIN_V
            pass
        pass
    
        if( block_face_grid_tag == GTAGS.MIN_V ) :
            if( hinge_dir_tag == DIR.U ) :
                hinge_grid_tag = GTAGS.CENTER_W
                ortho_grid_tag = GTAGS.MIN_U
                ortho_grid_select_tag = GTAGS.MAX_W
            pass
            if( hinge_dir_tag == DIR.W ) :
                hinge_grid_tag = GTAGS.CENTER_U
                ortho_grid_tag = GTAGS.MIN_W
                ortho_grid_select_tag = GTAGS.MIN_U
            pass
        pass
    
        if( block_face_grid_tag == GTAGS.MAX_V ) :
            if( hinge_dir_tag == DIR.U ) :
                hinge_grid_tag = GTAGS.CENTER_W
                ortho_grid_tag = GTAGS.MIN_U
                ortho_grid_select_tag = GTAGS.MIN_W
            pass
            if( hinge_dir_tag == DIR.W ) :
                hinge_grid_tag = GTAGS.CENTER_U
                ortho_grid_tag = GTAGS.MIN_W
                ortho_grid_select_tag = GTAGS.MAX_U
            pass
        pass
    
        if( block_face_grid_tag == GTAGS.MIN_W ) :
            if( hinge_dir_tag == DIR.U ) :
                hinge_grid_tag = GTAGS.CENTER_V
                ortho_grid_tag = GTAGS.MIN_U
                ortho_grid_select_tag = GTAGS.MIN_V
            pass
            if( hinge_dir_tag == DIR.V ) :
                hinge_grid_tag = GTAGS.CENTER_U
                ortho_grid_tag = GTAGS.MIN_V
                ortho_grid_select_tag = GTAGS.MAX_U
            pass
        pass
    
        if( block_face_grid_tag == GTAGS.MAX_W ) :
            if( hinge_dir_tag == DIR.U ) :
                hinge_grid_tag = GTAGS.CENTER_V
                ortho_grid_tag = +GTAGS.MIN_U
                ortho_grid_select_tag = GTAGS.MAX_V
            pass
            if( hinge_dir_tag == DIR.V ) :
                hinge_grid_tag = GTAGS.CENTER_U
                ortho_grid_tag = -GTAGS.MIN_V
                ortho_grid_select_tag = GTAGS.MIN_U
            pass
        pass
    
        return(hinge_grid_tag, ortho_grid_tag, ortho_grid_select_tag)
    pass


    #----------------------------------------------------------------------------------------------------
    # STATIC - GIVEN THE HINGE DIRECTION RETURN THE TAG THAT WILL RETURN THE START GRID ON THE HINGE AXIS
    @staticmethod
    def get_min_grid_tag(hinge_dir) : 
        min_grid_tag = -1
        if( hinge_dir == DIR.U ) :
            min_grid_tag = GTAGS.MIN_U
        pass
        if( hinge_dir == DIR.V ) :
            min_grid_tag = GTAGS.MIN_V
        pass
        if( hinge_dir == DIR.W ) :
            min_grid_tag = GTAGS.MIN_W
        pass

        return(min_grid_tag)
    pass

    #----------------------------------------------------------------------------------------------------
    # ASSOCIATE FACES BASED ON JOIN TYPE
    def sew(self) :
        #print '>>> TOP >>> ', here()

        ma = self.jface_a.block.mesh
        mb = self.jface_b.block.mesh

        #print 'NUMBER OF GRIDS IN MESH A =', len(ma.gl)
        #print 'NUMBER OF GRIDS IN MESH B =', len(mb.gl)

        #print 'A FACE TAG = ', self.jface_a.face_grid_tag
        #print 'B FACE TAG = ', self.jface_b.face_grid_tag

        self.jface_a.face_gl = ma.get_grid_list_from_tags(self.jface_a.face_grid_tag)
        self.jface_b.face_gl = mb.get_grid_list_from_tags(self.jface_b.face_grid_tag)
        
        #print '_'*80
        #print 'A FACE GRIDS:', self.jface_a.face_grid_tag, self.jface_a.name
        #c = -1
        #for ag in self.jface_a.face_gl :
        #  c += 1
        #  print c, ag
        #pass
    
        #print
        
        #print 'B FACE GRIDS:', self.jface_b.face_grid_tag, self.jface_b.name
        #c = -1
        #for bg in self.jface_b.face_gl :
        #  c += 1
        #  print c, bg
        #pass

        #print '_'*80



        ## FROM CALCULIX MANUAL...
        ##
        ##  The dependent surface is called the slave surface, the indepen-
        ## dent surface is the master surface. The user can freely decide which surface he
        ## takes as slave and which as master.
        ##
        ## Testing and the internet indicate that it is better if the slave surface is the denser of the two meshes
        ##


        
##    __ _______  ___________
##   / // /  _/ |/ / ___/ __/
##  / _  // //    / (_ / _/
## /_//_/___/_/|_/\___/___/
##

        if( self.type == JTAGS.HINGE ) :
            #print 'HINGE JOIN'

            a_ortho_grid_tag = -1
            b_ortho_grid_tag = -1
            a_ortho_grid_select_tag = -1
            b_ortho_grid_select_tag = -1
        

        
            a_hinge_grid_tag, a_ortho_grid_tag, a_ortho_grid_select_tag = join.get_face_hinge_line_tags(self.jface_a.face_grid_tag,
                                                                                                        self.jface_a.hinge_dir)
            self.jface_a.gl = ma.get_grid_list_from_tags(a_hinge_grid_tag,
                                                         self.jface_a.face_gl)
            
            self.jface_a.ortho_gl =  ma.get_grid_list_from_tags(a_ortho_grid_tag,
                                                                self.jface_a.face_gl)
                                                                

            b_hinge_grid_tag, b_ortho_grid_tag, b_ortho_grid_select_tag = join.get_face_hinge_line_tags(self.jface_b.face_grid_tag,
                                                                                                        self.jface_b.hinge_dir)
            self.jface_b.gl = mb.get_grid_list_from_tags(b_hinge_grid_tag,
                                                         self.jface_b.face_gl)
            
            self.jface_b.ortho_gl =  ma.get_grid_list_from_tags(b_ortho_grid_tag,
                                                                self.jface_b.face_gl)

            for g in self.jface_a.face_gl : 
                g.tags |= GTAGS.HINGED_FACE
            pass

            for g in self.jface_b.face_gl : 
                g.tags |= GTAGS.HINGED_FACE
            pass

            for g in self.jface_a.gl : 
                g.tags |= GTAGS.HINGED
            pass

            for g in self.jface_b.gl : 
                g.tags |= GTAGS.HINGED
            pass


            ## print 'NUMBER OF HINGE FACE GRIDS IN MESH A =', len(self.jface_a.face_gl)
            ## print 'NUMBER OF HINGE FACE GRIDS IN MESH B =', len(self.jface_b.face_gl)
            ## print 'NUMBER OF HINGE GRIDS IN MESH A =', len(self.jface_a.gl)
            ## print 'NUMBER OF HINGE GRIDS IN MESH B =', len(self.jface_b.gl)

            # SET THE SLAVE FACE TO BE THE ONE WITH THE MOST GRID POINTS ON THE HINGE LINE 
            if( len(self.jface_a.gl) >= len(self.jface_b.gl) ) :
                self.slave_tag = 'A'
                self.slave = self.jface_a
                self.master = self.jface_b
            else :
                self.slave_tag = 'B'
                self.slave = self.jface_b
                self.master = self.jface_a
            pass

        pass  # END JTAG.HINGE

##    ___  __  ______________
##   / _ )/ / / /_  __/_  __/
##  / _  / /_/ / / /   / /
## /____/\____/ /_/   /_/


        if( self.type == JTAGS.BUTT ) :
            #print 'BUTT JOIN'
            
            for g in self.jface_a.face_gl : 
                g.tags |= GTAGS.BUTT_FACE
            pass

            for g in self.jface_b.face_gl : 
                g.tags |= GTAGS.BUTT_FACE
            pass
        
            self.jface_a.get_face_lines()
            self.jface_b.get_face_lines()

            #print 'jface_a.face_gl LENGTH =', len(self.jface_a.face_gl)
            #print 'jface_b.face_gl LENGTH =', len(self.jface_b.face_gl)
            
            # SET THE SLAVE FACE TO BE THE ONE WITH THE MOST GRID POINTS ON THE FACE
            if( len(self.jface_a.face_gl) >= len(self.jface_b.face_gl) ) :
                self.slave_tag = 'A'
                self.slave = self.jface_a
                self.master = self.jface_b
            else :
                self.slave_tag = 'B'
                self.slave = self.jface_b
                self.master = self.jface_a
            pass
            
        pass # END JTAGS.BUTT

##    ___  ________________
##   / _ \/  _/ ___/  _/ _ \
##  / , _// // (_ // // // /
## /_/|_/___/\___/___/____/


        if( self.type == JTAGS.RIGID ) :
            #print 'RIGID JOIN'

            for g in self.jface_a.face_gl : 
                g.tags |= GTAGS.RIGID_FACE
            pass

            for g in self.jface_b.face_gl : 
                g.tags |= GTAGS.RIGID_FACE
            pass

            self.jface_a.get_face_lines()
            self.jface_b.get_face_lines()


            # SET THE SLAVE FACE TO BE THE ONE WITH THE MOST GRID POINTS ON THE FACE
            if( len(self.jface_a.face_gl) >= len(self.jface_b.face_gl) ) :
                self.slave_tag = 'A'
                self.slave = self.jface_a
                self.master = self.jface_b
            else :
                self.slave_tag = 'B'
                self.slave = self.jface_b
                self.master = self.jface_a
            pass
                  
        pass # END JTAGS.RIGID



    
        #print 'SLAVE GRIDS...'
        self.slave.face_el.clear()
        for sg in self.slave.face_gl :
            sg.tags |= GTAGS.JOIN_SLAVE
            #print sg
            eel = list(sg.el)
            for e in eel :
                #print '   USES ELEMENT ', e
                self.slave.face_el.add(e)
            pass
        pass

        #print 'MASTER GRIDS...'
        self.master.face_el.clear()
        for mg in self.master.face_gl :
            mg.tags |= GTAGS.JOIN_MASTER
            #print mg
            eel = list(mg.el)
            for e in eel :
                #print '   USES ELEMENT ', e
                self.master.face_el.add(e)
            pass
        pass


        #print 'SLAVE elements...'
        #tel = list(self.slave.face_el)
        #tel.sort(  key = lambda e : e.id )
        #for e in tel :
        #    print 'SLAVE ELEMENT # ', e.id
        #pass


        #print 'MASTER elements...'
        #tel = list(self.master.face_el)
        #tel.sort(  key = lambda e : e.id )
        #for e in tel :
        #    print 'MASTER ELEMENT # ', e.id
        #pass

        
    pass

pass

##########################################################################################################

# MODELS CAN HAVE MORE THAN ONE MESHES
class model(object) :

    def __init__ (self) :
        self.blocks = []
        self.meshes = []
        self.joins = []

        self.materials = []
        self.initial_temp = constants.DEFAULT_INITIAL_TEMP;
        self.reference_temp = constants.DEFAULT_REFERENCE_TEMP; # SEE CALCULIX *EXPANSION CARD
        self.NPSS_materials = {};
        self.results_elements = [] # THIS IS A LIST OF ELEMENT CREATED WHEN THE RESULTS ARE READ IF THEY ARE NOT FOUND IN THE MESHES

        self.safety_factor_yield = 1.0
        self.safety_factor_ultimate = 1.0
        self.header_text_lines = []
        self.do_nlgeom = True
        self.max_results = {}
        self.min_results = {}
        self.avg_results = {}
        self.eid_max = 0
        self.gid_max = 0        
    pass

    #----------------------------------------------------------------------------------------------------

    def __del__ (self) :
        del self.meshes
        del self.joins
        del self.materials
        del self.NPSS_materials
        del results_elements
        del self.header_text_lines
        del self.max_results
        del self.min_results
        del self.avg_results
    pass

    #----------------------------------------------------------------------------------------------------

    def plot_blocks(self) :
            
        fig = plt.figure()
        ax = fig.gca(projection='3d') # RETURNS CURRENT AXIS - CREATING ONE IF NECESSARY

        maxs = np.array([-1.0 * constants.BIG_REAL, -1.0 * constants.BIG_REAL, -1.0 * constants.BIG_REAL])
        mins = np.array([ constants.BIG_REAL, constants.BIG_REAL, constants.BIG_REAL])
        for b in self.blocks :
            #print ' LOOKING AT BLOCK :', b.name
            for p in b.pl.itervalues() :
                #print ' p = ', p
                for i in range(3) :
                    #print 'p.v['+str(i)+'] = ', p.v[i]
                    maxs[i] = max(p.v[i], maxs[i])
                    mins[i] = min(p.v[i], mins[i])
                pass
            pass
            b.plot_block_on_ax(ax)
        pass
       

        #print 'MAXS =', maxs
        #print 'MINS =', mins
        
        ax.legend()

        lo = min(mins) - 0.1 * abs(min(mins))
        hi = max(maxs) + 0.1 * abs(max(maxs))
        ax.set_xlim3d(lo, hi)
        ax.set_ylim3d(lo, hi)
        ax.set_zlim3d(lo, hi)

        ## vlo = mins - 0.1 * abs(mins)
        ## vhi = maxs + 0.1 * abs(maxs)
        ## ax.set_xlim3d(vlo[0], vhi[0])
        ## ax.set_ylim3d(vlo[1], vhi[1])
        ## ax.set_zlim3d(vlo[2], vhi[2])

        plt.show()
    pass

    #----------------------------------------------------------------------------------------------------
    # http://stackoverflow.com/questions/2709800/how-to-pickle-yourself

    def save(self, fname) :
        fp = open(fname, 'wb')
        #pickle.dump(self.__dict__, fp, 2)  # USING NEWER PROTOCOL
        pickle.dump(self.__dict__, fp, 0)  # OLDER ASCII BASED PROTOCOL
        fp.close()
    pass

    #----------------------------------------------------------------------------------------------------

    def read(self, fname):
        fp = open(fname,'rb')
        tmp_dict = pickle.load(fp)
        fp.close()          
        self.__dict__.update(tmp_dict)
    pass
    
    #----------------------------------------------------------------------------------------------------

    def set_initial_temp(self, temp) :
        self.initial_temp = temp
    pass

    #----------------------------------------------------------------------------------------------------

    def set_reference_temp(self, temp) :
        self.reference_temp = temp
    pass

    #----------------------------------------------------------------------------------------------------

    def add_block(self, bname) :
        b = block(bname, self)
        b.mesh = mesh(bname + '_MESH')
        b.mesh.block = b # BACK POINTER FROM THE MESH TO THIS BLOCK
        b.mesh.model = self # BACK POINTER FROM THE MESH TO THIS MODEL
        self.blocks.append(b)
        return(b)
    pass

    #----------------------------------------------------------------------------------------------------

    def add_material(self, mname) :
        m = material(mname)
        self.materials.append(m)
        return(m)
    pass

    #----------------------------------------------------------------------------------------------------

    def add_hinge_join(self, jname, ba, ba_face_grid_tag, ba_hinge_dir, bb, bb_face_grid_tag, bb_hinge_dir) :
        j = join(jname, ba, ba_face_grid_tag, bb, bb_face_grid_tag, JTAGS.HINGE, self)
        j.set_a_hinge_direction(ba_hinge_dir)
        j.set_b_hinge_direction(bb_hinge_dir)
        self.joins.append(j)
        return(j)
    pass

    #----------------------------------------------------------------------------------------------------

    def add_butt_join(self, jname, ba, ba_face_grid_tag, bb, bb_face_grid_tag ) :
        j = join(jname, ba, ba_face_grid_tag, bb, bb_face_grid_tag, JTAGS.BUTT, self)
        j.set_a_hinge_direction(None)
        j.set_b_hinge_direction(None)
        self.joins.append(j)
        return(j)
    pass

    #----------------------------------------------------------------------------------------------------

    def add_rigid_join(self, jname, ba, ba_face_grid_tag, bb, bb_face_grid_tag ) :
        j = join(jname, ba, ba_face_grid_tag, bb, bb_face_grid_tag, JTAGS.RIGID, self)
        j.set_a_hinge_direction(None)
        j.set_b_hinge_direction(None)
        self.joins.append(j)
        return(j)
    pass
    
    #----------------------------------------------------------------------------------------------------

    def generate_mesh(self) :
        for b in self.blocks :
            b.mesh.generate()
            pass
        pass
    pass

    #----------------------------------------------------------------------------------------------------

    def add_NPSS_material(self, mname) :
        m = None 
        if( len(self.NPSS_materials) == 0 ) :
            self.read_NPSS_materials()
        pass

        if( mname in self.NPSS_materials ) :
            #m = material(mname)
            m = copy.deepcopy(self.NPSS_materials[mname])
            # ASSIGN BACK POINTERS CORRECTLY IN COPIED MATERIAL
            for mp in m.mat_props :
                mp.material = m
            pass        
            self.materials.append(m)
        pass
        return(m)
    pass

    #----------------------------------------------------------------------------------------------------
    # http://www.faqs.org/docs/diveintopython/kgp_attributes.html
    
    def read_NPSS_materials(self) :


        s = os.path.dirname(__file__)
        ss = os.path.join(s, constants.MATS_FROM_NPSS)
        
        fp = open(ss)
        dom = minidom.parse(fp)
        fp.close()


        NPSS_required_fields = ['E', 'alpha', 'poissonRatio', 'rho']

        
        bad_material = False
        mat = None
        mp = None
        
        emats = dom.getElementsByTagName("MATERIAL")
        for em in emats :
            mat_name = None
            mat = None
            bad_material = False
            #print 'FOUND A MATERIAL'
            for k in em.attributes.keys() :
                aa = em.attributes[k]
                #print 'ATTRIBUTE = ', aa.name, ' VALUE =', aa.value
                if( aa.name == "name" ) :
                    mat_name = aa.value
                    #print 'GOT NAME ', mat_name
                    break
                pass
            pass
            if( mat_name is None ) :  # DIDN'T FIND A NAME FOR THE MATERIAL - SKIP IT
                continue
            pass
            mat = material(mat_name)
            # GET ALL THE TEMPERATURES FOR THE CURRENT MATERIAL
            etemps = em.getElementsByTagName("TEMPERATURE")
            for et in etemps :
                mat_temp = None
                for k in et.attributes.keys() :
                    aa = et.attributes[k]
                    #print ' ATTRIBUTE = ', aa.name, ' VALUE =', aa.value
                    if( aa.name == "temperature" ) :
                        mat_temp = aa.value
                        #print '   GOT TEMPERATURE ', mat_temp
                        break
                    pass
                pass

                if( mat_temp is None ) :  # DIDN'T FIND A TEMPERATURE VALUE - SKIP IT
                   continue
                pass

                mp = mat_prop()
                mp.temperature = float(mat_temp)
                
                # http://www.faqs.org/docs/diveintopython/kgp_child.html
                # http://www.faqs.org/docs/diveintopython/kgp_parse.html
                count_cards = dict.fromkeys(NPSS_required_fields, 0)

                eprops = [e for e in et.childNodes if e.nodeType == e.ELEMENT_NODE]
                for ep in eprops :
                    tag = ep.tagName
                    data = ep.firstChild.data
                    #print 'EP name =', ep.nodeName, ' TAG=', ep.tagName, ' DATA=|' + ep.firstChild.data + '|'
                    if( tag not in count_cards ) :
                        count_cards[tag] = 0
                    pass
                    if( data != "NG" ) :
                        count_cards[tag] = count_cards[tag] + 1
                        mp.props[tag] = float(data)
                    pass
                pass
                
                nn = 0
                for f in NPSS_required_fields :
                    nn = nn + count_cards[f]
                pass
                if( nn != len(NPSS_required_fields) ) :
                    bad_material = True
                    for p in mp.props :
                        del p
                    pass
                    mp.props.clear()
                    del mp
                    mp = None
                    break
                pass

                temp = float(mat_temp)
                ym = mp.props['E']
                tec = mp.props['alpha']
                pr = mp.props['poissonRatio']
                dens = mp.props['rho']
                sy = -1.0
                if( 'sy' in mp.props ) :
                    sy = mp.props['sy']
                pass
                su = -1.0
                if( 'su' in mp.props ) :
                    su = mp.props['su']
                pass
            
                mpp = mat.add_props_at_temp(temp, dens, ym, pr, tec, sy, su)

                # COPY OVER THE NPSS DATA FIELDS BUT MAKE SURE YOU DON'T OVERWITE ANY OF THE PROPERTIES USED HERE
                for k, v in mp.props.items() :
                    if( k not in mpp.props ) :
                        mpp.props[k] = v
                    pass
                pass


            pass # END TEMPERATURE LOOP
        
            if( bad_material ) :
                #print 'FOUND A BAD NPSS MATERIAL:', mat_name

                for pp in mat.mat_props :
                    del pp
                pass
                mat.mat_props = []
 
                del mat
                mat = None
                continue
            pass
        
            self.NPSS_materials[mat_name] = mat
            
        pass # END MATERIAL LOOP

        ## print 'NPSS MATERIALS RECOGNIZED:'
        ## cnt = 0
        ## for k, m  in self.NPSS_materials.items() :
        ##     cnt += 1
        ##     print '===================================='
        ##     print cnt, ") ", k
        ##     for mp in m.mat_props :
        ##         print '------------------------------ TEMP =', mp.temperature
        ##         keys = sorted(mp.props.keys())
        ##         for kk in keys :
        ##             print '    > ', kk, " = ", mp.props[kk]
        ##         pass
        ##     pass
        ## pass

            
    pass


    #----------------------------------------------------------------------------------------------------

    def read_stress_results_from_dot_frd(self, frd_fname) :

        frd_disp_begin_re = re.compile( r'^[ ]+[-][4][ ]+DISP[ ]+' )
        frd_stress_begin_re = re.compile( r'^[ ]+-4[ ]+STRESS[ ]+' )
        frd_field_tag_re = re.compile( r'^[ ]+[-][5][ ]+([A-Za-z0-9]+)[ ]+' )
        frd_data_line_start_re = re.compile( r'^[ ]+-1[ ]+')
        frd_data_line_next_re = re.compile( r'^[ ]+-2[ ]+')
        frd_data_line_end_re = re.compile( r'^[ ]+-3[ ]*$')

        fp = open(frd_fname, 'r')
        lines = [ line.replace('\r\n', '') for line in fp ]
        fp.close()

        el = []
        for b in self.blocks :
            el = el + b.mesh.el.values()
        pass
    
        gl = []
        for b in self.blocks :
            gl = gl + b.mesh.gl.values()
        pass
        gl = filter( ( lambda g: g is not None ), gl)


        # RIGHT NOW WE ARE ASSUMING WE ARE READING THE RESULTS OF THE MESH THAT IS IN MEMORY
        #  IN THE FUTURE WE COULD (SHOULD?) BE ABLE TO READ THE FULL MESH IN FROM THE FRD FILE
        #    RIGHT NOW WE WILL ASSUME THAT IS REDUNDANT.

        # NOTE: FOR NLGEOM ANLAYSES WITH MUTLIPLE TIME INCREMENTS
        #  THE DISPLACMENTS AND STRESSES WILL READ FOR EACH TIME INCREMENT
        #  OVERWRITING ANY VALUES IN MEMORY FROM EARLIER INCREMENTS
        #  THE LAST SET OF DISPLACEMENTS AND STRESSES WILL BE USED
        #  TO FIND THE MAX AND MIN VALUES.

        lit = iter(lines)

        fw = 12 # FIXED FIELD WIDTH FOR DISPLACEMENTS AND STRESSES

        for i in range( len(lines) ) :
            try:
                l = lit.next()
                
                # DISPLACEMENTS
                dt = frd_disp_begin_re.search(l)
                if( dt ) :
                    #print '************** FOUND DISPLACEMENTS'
                    fields = []
                    fields.append('ID')
                    while( not frd_data_line_start_re.search(l) ) :
                        ft = frd_field_tag_re.search(l)
                        if( ft is not None ) :
                            fields.append(ft.group(1))
                        pass
                        l = lit.next()
                    pass
                    #print 'FIELDS = ', fields
                    fids = enum(fields)
                    while True :
                        vals = []
                        l = re.sub(r'^ -1 ', '    ', l)
                        for ff in range(4) :
                            fs = ff * fw + 1
                            fe = fs + fw
                            vals.append( float( l[fs:fe].strip() ) )
                            if( ff == 0 ) :
                                vals[ff] = int(vals[ff])
                            pass
                        pass
                        all = math.sqrt(vals[1]*vals[1] + vals[2]*vals[2] + vals[3]*vals[3])
                        vals.append(all)
                        

                        res = {}
                        #print 'DISP ',
                        for ff in fields :
                            iid = eval("fids." + ff)
                            #print vals[iid], ' ',
                            res[ff] = vals[iid]
                        pass
                        #print
                       
                        gegl = filter( ( lambda gg: gg.id == vals[fids.ID] ), gl)
                        if( len(gegl) == 0 ) :
                            print 'READING FRD RESULTS - DISPLACEMENTS: COULD NOT FIND GRID: ', id
                            print ' ---- SKIPPING GRID RESULTS'
                        pass
                    
                        if( len(gegl) > 1 ) :
                            print 'READING FRD RESULTS - DISPLACEMENTS: MORE THAN ONE GRID OBJECTS FOR GRID: ', id, ' >', len(gegl)
                            print ' ---- ASSIGNING GRID RESULTS TO ALL GRIDS FOUND'
                        pass

                        for gp in gegl :
                            for k, v in res.items() :
                                gp.results[k] = v;
                            pass
                        pass

                        l = lit.next()
                        if( frd_data_line_end_re.search(l) ) :
                            break
                        pass
                    pass
                pass

                # STRESS
                dt = frd_stress_begin_re.search(l)
                if( dt ) :
                    fields = []
                    fields.append('ID')
                    while( not frd_data_line_start_re.search(l) ) :
                        ft = frd_field_tag_re.search(l)
                        if( ft is not None ) :
                            fields.append(ft.group(1))
                        pass
                        l = lit.next()
                    pass
                    #print 'FIELDS = ', fields
                    fids = enum(fields)
                    while True :
                        vals = []
                        l = re.sub(r'^ -1 ', '    ', l)
                        for ff in range(7) :
                            fs = ff * fw + 1
                            fe = fs + fw
                            vals.append( float( l[fs:fe].strip() ) )
                            if( ff == 0 ) :
                                vals[ff] = int(vals[ff])
                            pass
                        pass

                        res = {}
                        #print 'STRESS ',
                        for ff in fields :
                            iid = eval("fids." + ff)
                            #print vals[iid], ' ',
                            res[ff] = vals[iid]
                        pass
                        #print
                        
                        dsxy = res['SXX'] - res['SYY']
                        dsyz = res['SYY'] - res['SZZ']
                        dsxz = res['SXX'] - res['SZZ']
                        dss = ( res['SYZ'] * res['SYZ']
                                + res['SZX'] * res['SZX']
                                + res['SXY'] * res['SXY'] )
                        
                        vm = math.sqrt( 1.0 / 2.0 * ( dsxy * dsxy + dsyz * dsyz + dsxz * dsxz + 6.0 * dss ) )
                        res['VM'] = vm
                        #if( vm < minvm ) :
                        #    print 'NEW MIN VM = ', vm, ' OLD MIN VM WAS = ', minvm
                        #    minvm = vm
                        #    print 'RES=', res
                        #pass

                        #print 'STRESS: LOOKING FOR GRID WITH GID = ', vals[fids.ID]
                        gegl = filter( ( lambda gg: gg.id == vals[fids.ID] ), gl)
                        #print 'FOUND GRIDS...', gegl
                        
                        if( len(gegl) == 0 ) :
                            print 'READING FRD RESULTS - STRESSES: COULD NOT FIND GRID OBJECT FOR GRID: ', id
                            print ' ---- SKIPPING GRID RESULTS'
                        pass
                    
                        if( len(gegl) > 1 ) :
                            print 'READING FRD RESULTS - STRESSES: MORE THAN ONE GRID OBJECTS FOR GRID: ', id, ' >', len(gegl)
                            print ' ---- ASSIGNING GRID RESULTS TO ALL GRIDS FOUND'
                        pass

                        for gp in gegl :
                            #print 'ASSIGNING RESULTS VALUES TO GRID:', gp.id 
                            for k, v in res.items() :
                                #print 'KEY=', k, '   VALUE=', v
                                gp.results[k] = v;
                            pass
                        pass

                        l = lit.next()
                        if( frd_data_line_end_re.search(l) ) :
                            break
                        pass
                    pass
                pass

                continue

            except StopIteration, e:
                #print 'STOP ITERATION'
                break
            except :
                print "Unexpected error while reading FRD file:", sys.exc_info()[0]
                raise
            pass
        pass



        # FOR THE GRIDS THE MAX, MIN AND AVG RESULTS ARE THE SAME AS THE RESULTS THAT WERE READ IN
        for g in gl :
            for k, v in g.results.items() :
                g.max_results[k] = v
                g.min_results[k] = v
                g.avg_results[k] = v
            pass
        pass

        # CALC ELEMENT RESULTS : MEMBER GRID MAX/MIN/AVG
        for e in el :
            g = e.avg_grid
            g.results = {}
            g.max_results = {}
            g.min_results = {}
            g.avg_results = {}

            nexg = len(e.gl)
            if( nexg == 0 ) : continue
            for eg in e.gl :
                
                for k, v in eg.avg_results.items() :
                    # AVGERAGE RESULTS
                    if( k not in g.avg_results ) :
                        g.avg_results[k] = v / float(nexg)
                    else :
                        g.avg_results[k] += ( v / float(nexg) ) 
                    pass
                pass

                for k, v in eg.max_results.items() :
                    # MAX RESULTS
                    if( k not in g.max_results ) :
                        g.max_results[k] = v
                    pass
                    g.max_results[k]= max(g.max_results[k], v)
                pass

                for k, v in eg.min_results.items() :
                    # MIN RESULTS
                    if( k not in g.min_results ) :
                        g.min_results[k] = v
                    pass
                    g.min_results[k]= min(g.min_results[k], v)
                pass
            pass

        pass
    
        
        # CALCULATE MOS's IF WE HAVE A MESH WITH MATERIAL PROPS - BASED ON VONMISES
        # http://en.wikipedia.org/wiki/Factor_of_safety
        for e in el :
            if( e.mesh is not None ) :
                temp = e.avg_grid.get_param('TEMP', self.initial_temp)
                vf = e.avg_grid.get_param('VOL_FRAC', 1.0)
                tmatprop = e.mesh.material.get_mat_props_at_temp_with_volume_fraction_scaling(temp, vf)
                #print 'BEFORE calc_mos AVG GRID'
                self.calc_mos(e.avg_grid, tmatprop)
                for eg in e.gl :
                    #print 'BEFORE calc_mos element grid = ', eg.id
                    self.calc_mos(eg, tmatprop)
                pass
                    
            pass

        pass


        # FIND OVERALL MAX AND MINS:
        #   NOW LOOK AT ALL THE GRIDS AND FIND THE MAX AND MINS
        tgl = gl
        for m in self.blocks :
            for e in b.mesh.el.values() :
                tgl.append(e.avg_grid)
            pass
        pass
    

        self.max_results = {}
        self.min_results = {}
        self.avg_results = {}

        for g in tgl :
            
            for k, v in g.avg_results.items() :        
                # AVGERAGE RESULTS - THIS IS PRETTY WORTHLESS
                if( k not in self.avg_results ) :
                    self.avg_results[k] = v / float(len(tgl))
                else :
                    self.avg_results[k] += ( v / float(len(tgl)) )
                pass
            pass
        
            for k, v in g.max_results.items() :
                #if( k == 'SY' ) :
                #    print 'looking at a SY RESULT FOR GRID #', g.id, ' SY=', v
                #pass
                # MAX RESULTS
                if( k not in self.max_results ) :
                    self.max_results[k] = v
                    self.max_results[k + '_ID'] = g.id
                    
                    #if( k == 'SY' ) :
                    #    print 'FIRST MAX', self.max_results[k], ' ID=', self.max_results[k + '_ID']
                    #pass
                
                pass
                if( self.max_results[k] < v ) :
                    self.max_results[k] = v
                    self.max_results[k + '_ID'] = g.id                    
                    #if( k == 'SY' ) :
                    #    print 'NEW MAX', self.max_results[k], ' ID=', self.max_results[k + '_ID']
                    #pass
                pass
            pass
        
            for k, v in g.min_results.items() :
                # MIN RESULTS
                if( k not in self.min_results ) :
                    self.min_results[k] = v
                    self.min_results[k + '_ID'] = g.id
                pass
                if( self.min_results[k] > v ) :
                    self.min_results[k] = v
                    self.min_results[k + '_ID'] = g.id
                pass
            pass
        pass


    pass

    #----------------------------------------------------------------------------------------------------

    def calc_mos(self, g, matp) :
        #print matp.props
        #print 'VM MIN= ', g.min_results['VM'], '  AVG= ', g.avg_results['VM'], '  MAX= ', g.max_results['VM']
        # ULTIMATE MOS

        # GRIDS ARE SHARED BY ELEMENTS THAT MAY HAVE DIFFERENT MATERIAL PROPERTIES
        # STORE AT THE GRID THE WORST CASE MOSU AND MOSY FOR ALL ELEMENTS THAT REFERENCE A GIVEN GRID

        if( matp.props['SU'] > 0.0 ) :
            
            # AVERAGE - REALLY THE MINIMUM OF THE AVERAGE - Not much use really
            if(  g.avg_results['VM'] > constants.TOL ) :
                mosu = matp.props['SU'] / ( g.avg_results['VM'] * self.safety_factor_ultimate ) - 1.0
            else :
                mosu = constants.BIG_REAL
            pass
            update = True
            #print 'AVG MOSU = ', mosu, 'VM = ', g.avg_results['VM']
            # KEEP THE SMALLEST MOS FOR EACH GRID - ITS THE ONE THAT CONTROLS
            if( 'MOSU' in g.avg_results ) :
                if( mosu > g.avg_results['MOSU'] ) :
                    update = False
                pass
            pass
            if( update ) :
                g.avg_results['MOSU'] = mosu
                g.avg_results['SU'] = matp.props['SU']
                g.avg_results['SFU'] = self.safety_factor_ultimate
            pass

            # MAX STRESS WILL GIVE THE MINIMUM MOS
            if( g.min_results['VM']  > constants.TOL ) :
                mosu_mins = matp.props['SU'] / ( g.min_results['VM'] * self.safety_factor_ultimate ) - 1.0
            else :
                mosu_mins = constants.BIG_REAL
            pass
        
            if( g.max_results['VM']  > constants.TOL ) :
                mosu_maxs = matp.props['SU'] / ( g.max_results['VM'] * self.safety_factor_ultimate ) - 1.0
            else :
                mosu_maxs = constants.BIG_REAL
            pass
            
            # MIN
            mosu = min(mosu_mins, mosu_maxs)
            update = True
            #print 'MIN MOSU = ', mosu
            # KEEP THE SMALLEST MOS FOR EACH GRID
            if( 'MOSU' in g.min_results ) :
                if( mosu > g.min_results['MOSU'] ) :
                    update = False
                pass
            pass
            if( update ) :
                g.min_results['MOSU'] = mosu
                g.min_results['SU'] = matp.props['SU']
                g.min_results['SFU'] = self.safety_factor_ultimate
            pass


            # MAX
            mosu = max(mosu_mins, mosu_maxs)
            update = True
            #print 'MAX MOSU = ', mosu
            # KEEP THE LARGEST MOS FOR EACH GRID 
            if( 'MOSU' in g.max_results ) :
                if( mosu < g.max_results['MOSU'] ) :
                    update = False
                pass
            pass
            if( update ) :
                g.max_results['MOSU'] = mosu
                g.max_results['SU'] = matp.props['SU']
                g.max_results['SFU'] = self.safety_factor_ultimate
            pass

        # YEILD MOS
        if( matp.props['SY'] > 0.0 ) :
            
            # AVERAGE - REALLY THE MINIMUM OF THE AVERAGE
            if( g.avg_results['VM']  > constants.TOL ) :
                mosy = matp.props['SY'] / ( g.avg_results['VM'] * self.safety_factor_yield ) - 1.0
            else :
                mosy = constants.BIG_REAL
            pass
        
            update = True
            # KEEP THE SMALLEST MOS FOR EACH GRID - ITS THE ONE THAT CONTROLS
            if( 'MOSY' in g.avg_results ) :
                if( mosy > g.avg_results['MOSY'] ) :
                    update = False
                pass
            pass
            if( update ) :
                g.avg_results['MOSY'] = mosy
                g.avg_results['SY'] = matp.props['SY']
                g.avg_results['SFY'] = self.safety_factor_yield
            pass

            # MAX STRESS WILL GIVE THE MINIMUM MOS
            if( g.min_results['VM']  > constants.TOL ) :
                mosy_mins = matp.props['SY'] / ( g.min_results['VM'] * self.safety_factor_yield ) - 1.0
            else :
                mosy_mins = constants.BIG_REAL
            pass
        
            if( g.max_results['VM']  > constants.TOL ) :
                mosy_maxs = matp.props['SY'] / ( g.max_results['VM'] * self.safety_factor_yield ) - 1.0
            else :
                mosy_maxs = constants.BIG_REAL
            pass
        
            # MIN
            mosy = min(mosy_mins, mosy_maxs)
            update = True
            # KEEP THE SMALLEST MOS FOR EACH GRID
            if( 'MOSY' in g.min_results ) :
                if( mosy > g.min_results['MOSY'] ) :
                    update = False
                pass
            pass
            if( update ) :
                g.min_results['MOSY'] = mosy
                g.min_results['SY'] = matp.props['SY']
                g.min_results['SFY'] = self.safety_factor_yield
            pass

            # MAX
            mosy = max(mosy_mins, mosy_maxs)
            update = True
            # KEEP THE LARGEST MOS FOR EACH GRID 
            if( 'MOSY' in g.max_results ) :
                if( mosy < g.max_results['MOSY'] ) :
                    update = False
                pass
            pass
            if( update ) :
                g.max_results['MOSY'] = mosy
                g.max_results['SY'] = matp.props['SY']
                g.max_results['SFY'] = self.safety_factor_yield
            pass

        pass
        #print g.max_results
    pass
        

    #----------------------------------------------------------------------------------------------------

    def stitch(self) :
                          
        for j in self.joins :
            #print 'JOIN : ', j
            j.sew()
        pass
    pass

    #----------------------------------------------------------------------------------------------------

    def save_abq(self, fname) :
        #print '>>> TOP >>> ', here()

        face_tags = [ FTAGS.MAX_U, FTAGS.MIN_U, FTAGS.MAX_V, FTAGS.MIN_V, FTAGS.MAX_W, FTAGS.MIN_W ]
        face_tags_exploded = [ FTAGS.MAX_W, FTAGS.MIN_W ]

        self.stitch()

        b2b_brace_re = re.compile( r'[}][ ]*[{]', re.IGNORECASE )

        # LIST OF ALL GRIDS
        ggl = []
        for b in self.blocks :
            m = b.mesh
            tgl = m.gl.values()
            tgl = filter( ( lambda g: g is not None ), tgl)
            tgl.sort(key=lambda g : g.id)
            ggl = ggl + tgl
        pass
        ggl.sort(key=lambda g : g.id)

        #for g in ggl :
        #    print 'GRID:', g
        #pass
        

        # LIST OF ALL ELEMENTS
        gel = []
        for b in self.blocks :
            m = b.mesh
            tel = m.el.values()
            tel.sort(key=lambda e : e.id)
            gel = gel + tel
        pass
        gel.sort(key=lambda e : e.id)
        #for e in gel :
        #    print 'ELEM:', e
        #pass



        for m in self.materials :
            m.mat_props.sort(key = lambda prp: prp.temperature)
            m.lo_temp_defined = m.mat_props[0].temperature
            m.hi_temp_defined = m.mat_props[-1].temperature
        pass

        for b in self.blocks :
            m = b.mesh
            m.calc_element_parameter_averages()
        pass

        now_here_time = now()

        fp = open(fname, 'w+')
        fp.write('*HEADING' + '\n')
        fp.write('MODEL DYNAMICALLY GENERATED ' + fname + '\n')
        fp.write('** THIS MESH WAS GENERATED ON ' + now_here_time + '\n')
        fp.write('*' * 80 + '\n')

        for l in self.header_text_lines :
            fp.write('** ' + l + '\n')
        pass
    
        fp.write('*' * 80 + '\n')

        nn = []
        for b in self.blocks :
            m = b.mesh
            nn.append(m.nset_name)      
            fp.write('** <BEGIN> ' + '~-~' * 15 + ' ' + m.nset_name  +' NODES\n')
            fp.write('*NODE, NSET=' + m.nset_name + '\n')
            m.save_abq_grids(fp)
            fp.write('** <END> ' + '~-~' * 15 + ' ' + m.nset_name  +' NODES\n')
        pass

        fp.write('*NSET, NSET=NALL' + '\n')
        
        tc = 0
        lc = 0
        for ns in nn :
            tc = tc + 1
            lc = lc + 1
            ns = ns.strip()
            fp.write(ns)
            if( ( lc == 13 ) or ( tc == len(nn) ) ) :
                lc = 0
                if( tc != len(nn) ) :
                    fp.write(', ')
                pass
                fp.write('\n')
            else :
                fp.write(', ')
            pass
        pass


        for b in self.blocks :
            m = b.mesh
            fp.write('** <BEGIN> ' + '~-~' * 15 + ' ' + m.eset_name  +' ELEMENTS\n')
            m.save_abq_elements(fp)
            fp.write('** <END> ' + '~-~' * 15 + ' ' + m.eset_name  +' ELEMENTS\n')
        pass


        fp.write('*ELSET, ELSET=EALL\n')
        ii = 0
        ti = 0
        for e in gel :
            ii += 1
            ti += 1
            fp.write(str(e.id))
            if( ii == 13) :
                fp.write('\n')
                ii = 0
            else :
                if( ti != len(gel) ) : 
                    fp.write(', ')
                pass
            pass
        pass
        if( ii != 0 ) :
            fp.write('\n')
        pass

        fp.write('*' * 80 + '\n')



        for b in self.blocks :
            m = b.mesh
            #print 'MESH ', m.name, ' USES MATERIAL ', m.material.name
            #print '   TLO=', m.temp_lo, '  THI=', m.temp_hi
            lot = m.temp_lo * 0.95
            hit = m.temp_hi * 1.05
            if( abs(hit - lot) < 1.0 ) :
                lot = lot - 1.0
                hit = hit + 1.0
            pass
            dt = (hit - lot) / m.material.number_of_output_divisions
            tel = m.el.values()
            tel.sort(key=lambda ee : ee.id)
            for e in tel :
                vf = e.avg_grid.get_param('VOL_FRAC', 1.0)
                
                #print 'ELEMENT VOLUME FRACTION = ', vf
                #print 'ELE AVG GRID PARAMS = ', e.avg_grid.params
                
                ml = []
                for it in range(m.material.number_of_output_divisions + 1) :
                    t = lot + float(it) * dt
                    tmatprop = m.material.get_mat_props_at_temp_with_volume_fraction_scaling(t, vf)
                    ml.append(tmatprop)
                pass
                mat_name = m.material.name + '_' + e.name
                fp.write('*MATERIAL, NAME=' + mat_name + '\n')
                fp.write('**  ELSET VOL FRAC = ' + str(vf) + '\n')
                fp.write('*DENSITY' + '\n')
                f = '{0}, {1}'
                for mat in ml :
                    fp.write(f.format(mat.props['DENSITY'], mat.temperature) + '\n')
                pass
                fp.write('*ELASTIC, TYPE=ENGINEERING CONSTANTS' + '\n')
                
                f1 = b2b_brace_re.sub(r'}, {',  ' {}' * 8).strip()
                f2 = b2b_brace_re.sub(r'}, {',  ' {}' * 2).strip()
                for mat in ml :
                    s = f1.format(mat.props['YM11'],
                                  mat.props['YM22'],
                                  mat.props['YM33'],
                                  mat.props['PR12'],
                                  mat.props['PR13'],
                                  mat.props['PR23'],
                                  mat.props['SM12'],
                                  mat.props['SM13'])
                    fp.write(s + '\n')
                    s = f2.format(mat.props['SM23'],
                                  mat.temperature)
                    fp.write(s + '\n')
                pass
                fp.write('*EXPANSION, TYPE=ORTHO, ZERO=' + str(self.reference_temp) + '\n')
                f = b2b_brace_re.sub(r'}, {',  ' {}' * 4).strip()
                for mat in ml :
                    s = f.format(mat.props['TEC11'],
                                 mat.props['TEC22'],
                                 mat.props['TEC33'],
                                 mat.temperature)
                    fp.write( s + '\n' )
                pass

                s = '*SOLID SECTION, ELSET=' + e.single_element_set_name + ', MATERIAL=' + mat_name 
                fp.write(s + '\n')
            pass
            fp.write('*' * 80 + '\n')
        pass



        fp.write('*' * 80 + '\n')
        feq = '{}, S{}'
        
        for j in self.joins :
            
            jfs = j.slave
            jfm = j.master            

            if( j.type == JTAGS.HINGE ) :
                slave_surf_name = 'S_' + jfs.name + '_LINE'
                master_surf_name = 'M_' + jfm.name + '_SURF'

                # WRITE SLAVE NODE SET
                fp.write('** <BEGIN> ' + '~-~' * 5 + ' ' + jfs.name  +' HINGE LINE NODE SET\n')
                self.write_nset(fp, jfs.name, jfs.gl)
                fp.write('** <END>   ' + '~-~' * 5 + ' ' + jfs.name  +' HINGE LINE NODE SET\n')

                # WRITE SLAVE LINE SURFACE
                fp.write('** <BEGIN> ' + '~-~' * 5 + ' ' + slave_surf_name  +' HINGE SLAVE SURFACE (LINE)\n')
                fp.write('*SURFACE, NAME=' + slave_surf_name + ', TYPE=NODE\n')
                fp.write(jfs.name + '\n')
                fp.write('** <END>   ' + '~-~' * 5 + ' ' + slave_surf_name +' HINGE SLAVE SURFACE (LINE)\n')

                # WRITE MASTER SURFACE
                fp.write('** <BEGIN> ' + '~-~' * 5 + ' ' + master_surf_name  +' HINGE MASTER SURFACE\n')
                fp.write('*SURFACE, NAME=' + master_surf_name + ', TYPE=ELEMENT\n')
                
                tel = list(jfm.face_el)
                tel.sort(  key = lambda e : e.id )
                face_num = constants.FACE_NUM(jfm.face_grid_tag)
                #face_num = int( ( math.log(jfm.face_grid_tag) / math.log(2.0) ) + 1 )
                for me in tel :
                    fp.write(feq.format(me.id, face_num) + '\n')
                pass
                fp.write('** <END>   ' + '~-~' * 5 + ' ' + master_surf_name  +' HINGE MASTER SURFACE\n')
                
                fp.write('*TIE \n')
                fp.write(slave_surf_name + ', ' + master_surf_name + '\n')
                fp.write('*' * 80 + '\n')
               
            pass


            if( j.type == JTAGS.RIGID ) :

                for jf in j.jfaces :
                    fp.write('** <BEGIN> ' + '~-~' * 15 + ' ' + jf.name  +' JOIN NODE SET\n')
                    self.write_nset(fp, jf.name, jf.face_gl)
                    fp.write('** <END>   ' + '~-~' * 15 + ' ' + jf.name  +' JOIN NODE SET\n')
                pass
                fp.write('*NSET, NSET=' + j.name + '\n')
                s = ', '.join([ jf.name for jf in j.jfaces ])
                fp.write(s + '\n')

                fp.write('*RIGID BODY, NSET=' + j.name +'\n')
                fp.write('*' * 80 + '\n')
            pass

            
            if( j.type == JTAGS.BUTT ) :

                fp.write('** <BEGIN> ' + '~-~' * 5 + ' ' + jfs.name  +' JOIN SLAVE SURFACE\n')
                fp.write('*SURFACE, NAME=' + jfs.name + 'S, TYPE=ELEMENT\n')
                
                # FROM CALCULIX MANUAL... face_num
                #  face 1: 1-2-3-4 -> FTAGS.MIN_W =  1 = GTAGS.MIN_W
                #  face 2: 5-8-7-6 -> FTAGS.MAX_W =  2 = GTAGS.MAX_W
                #  face 3: 1-5-6-2 -> FTAGS.MIN_V =  4 = GTAGS.MIN_V
                #  face 4: 2-6-7-3 -> FTAGS.MAX_U =  8 = GTAGS.MAX_U
                #  face 5: 3-7-8-4 -> FTAGS.MAX_V = 16 = GTAGS.MAX_V
                #  face 6: 4-8-5-1 -> FTAGS.MIN_U = 32 = GTAGS.MIN_U

                tel = list(jfs.face_el)
                tel.sort(  key = lambda e : e.id )
                face_num = constants.FACE_NUM(jfs.face_grid_tag)
                #face_num = int( ( math.log(jfs.face_grid_tag) / math.log(2.0) ) + 1 )
                for se in tel :
                    fp.write(feq.format(se.id, face_num) + '\n')
                pass

                fp.write('** <END>   ' + '~-~' * 5 + ' ' + jfs.name  +' JOIN SLAVE SURFACE\n')
            
                fp.write('** <BEGIN> ' + '~-~' * 5 + ' ' + jfm.name  +' JOIN MASTER SURFACE\n')
                fp.write('*SURFACE, NAME=' + jfm.name + 'S, TYPE=ELEMENT\n')
                
                # FROM CALCULIX MANUAL... face_num
                #  face 1: 1-2-3-4 -> FTAGS.MIN_W =  1 = GTAGS.MIN_W
                #  face 2: 5-8-7-6 -> FTAGS.MAX_W =  2 = GTAGS.MAX_W
                #  face 3: 1-5-6-2 -> FTAGS.MIN_V =  4 = GTAGS.MIN_V
                #  face 4: 2-6-7-3 -> FTAGS.MAX_U =  8 = GTAGS.MAX_U
                #  face 5: 3-7-8-4 -> FTAGS.MAX_V = 16 = GTAGS.MAX_V
                #  face 6: 4-8-5-1 -> FTAGS.MIN_U = 32 = GTAGS.MIN_U

                tel = list(jfm.face_el)
                tel.sort(  key = lambda e : e.id )
                face_num = constants.FACE_NUM(jfm.face_grid_tag)
                #face_num = int( ( math.log(jfm.face_grid_tag) / math.log(2.0) ) + 1 )
                for me in tel :
                    fp.write(feq.format(me.id, face_num) + '\n')
                pass
            
                fp.write('** <END>   ' + '~-~' * 5 + ' ' + jfm.name  +' JOIN MASTER SURFACE\n')


                fp.write('*TIE \n')
                fp.write(jfs.name + 'S, '+ jfm.name + 'S' + '\n')

                fp.write('*' * 80 + '\n')

            
            pass
        pass




        fp.write('*INITIAL CONDITIONS, TYPE=TEMPERATURE' + '\n')
        fp.write('NALL, ' + str(self.initial_temp) + '\n' )
        

        fp.write('*' * 80 + '\n')


        fp.write('*BOUNDARY' + '\n')
        for b in self.blocks :
            m = b.mesh
            tgl = m.gl.values()
            tgl = filter( ( lambda g: g is not None ), tgl)
            tgl.sort(key=lambda g : g.id)            
            for g in tgl :
                if( g.id <= 0 ) : continue
                #print 'gid = ', g.id, ' bc\'s =', g.bc
                bcl = sorted(list(g.bc))
                if( len(bcl) > 0 ) :
                    for dof in bcl :
                        fp.write(str(g.id) + ', ' + str(dof) + '\n')
                    pass
                pass
                #fp.write(str(g.id) + ', 11, 11, ' + str(g.temp) + '\n') 
            pass
        pass

 
        if( self.do_nlgeom ) :
            fp.write('*STEP, INC=100, NLGEOM' + '\n')
            fp.write('*STATIC' + '\n')
            fp.write('0.1, 1.0' + '\n')
        else :
            fp.write('*STEP' + '\n')
            fp.write('*STATIC' + '\n')
        pass
    

        fp.write('*DLOAD' + '\n')

        vl = [None] * 3 # WORKING VECTOR LIST


        #print 'WRITE PRESSURE'
        for b in self.blocks :
            m = b.mesh
            # IF THIS IS AN EXPLODED BLOCK ONLY PUT PRESSURES ON THE +/- W FACES
            if( b.defined_by_tag == BTAGS.EXPLODED ) :
                face_tags = face_tags_exploded
            pass
            tel = m.el.values()
            tel.sort(  key = lambda e : e.id )

            #for e in tel :
            #    print 'ELEM:', e
            #pass


            #print 'FACE_TAGS =', face_tags
            
            for ft in face_tags :
                #print 'FACETAG = ', ft
                face_num = constants.FACE_NUM(ft)
                #print 'LOOKING AT FACENUM = ', face_num

                ftel = filter( ( lambda ee : ee.tags & ft ), tel )
                #ftel.sort(  key = lambda e : e.id )

                # WE DO ASSUME THAT THE START AND END POINTS OF NEIGHBORING EDGES
                # OF A GIVE BLOCK FACE DEFINITION ARE IDENTICAL, SO WE ONLY REFRENCE THE
                # START POINT HERE.
                # THIS IS AN ARTIFACT OF HOW WE ALWAYS BUILD THE BLOCK FACES.
                bface = b.faces[ft]
                if( ft == FTAGS.MIN_U ) :
                    fnn = bface.edges[3].sp
                    fpn = bface.edges[0].sp
                    fnp = bface.edges[2].sp
                    fpp = bface.edges[1].sp
                    # V -> W
                    interp_dirs = [ DIR.V, DIR.W ]
                pass
                if( ft == FTAGS.MAX_U ) :
                    fnn = bface.edges[0].sp
                    fpn = bface.edges[3].sp
                    fnp = bface.edges[1].sp
                    fpp = bface.edges[2].sp
                    # V -> W
                    interp_dirs = [ DIR.V, DIR.W ]
                pass

            
                if( ft == FTAGS.MIN_V ) :
                    fnn = bface.edges[0].sp
                    fpn = bface.edges[3].sp
                    fnp = bface.edges[1].sp
                    fpp = bface.edges[2].sp
                    # U -> W
                    interp_dirs = [ DIR.U, DIR.W ]
                pass
                if( ft == FTAGS.MAX_V ) :
                    fnn = bface.edges[3].sp
                    fpn = bface.edges[0].sp
                    fnp = bface.edges[2].sp
                    fpp = bface.edges[1].sp
                    # U -> W
                    interp_dirs = [ DIR.U, DIR.W ]
                pass
            
                if( ft == FTAGS.MIN_W ) :
                    fnn = bface.edges[0].sp
                    fpn = bface.edges[1].sp
                    fnp = bface.edges[3].sp
                    fpp = bface.edges[2].sp
                    # U -> V
                    interp_dirs = [ DIR.U, DIR.V ]
                pass
                if( ft == FTAGS.MAX_W ) :
                    fnn = bface.edges[0].sp
                    fpn = bface.edges[3].sp
                    fnp = bface.edges[1].sp
                    fpp = bface.edges[2].sp
                    # U -> V
                    interp_dirs = [ DIR.U, DIR.V ]
                pass

                #print 'fnn pressure = ',  fnn.get_param('PRES', 0.0)
                #print 'fpn pressure = ',  fpn.get_param('PRES', 0.0)
                #print 'fnp pressure = ',  fnp.get_param('PRES', 0.0)
                #print 'fpp pressure = ',  fpp.get_param('PRES', 0.0)

                psum = 0.0
                psum += abs(fnn.get_param('PRES', 0.0))
                psum += abs(fpn.get_param('PRES', 0.0))
                psum += abs(fnp.get_param('PRES', 0.0))
                psum += abs(fpp.get_param('PRES', 0.0))
                psum = psum / 4.0
                if( abs(psum) < constants.TOL ) : continue

                #  TBD - WE MIGHT WANT TO IMPLIMENT SOMETHING LIKE THE FOLLOWING...
                #                
                # #http://stackoverflow.com/questions/5124126/python-scipy-interpolation-map-coordinates#
                # import numpy
                # from scipy import interpolate
                # x = numpy.array([0.0, 0.60, 1.0])
                # y = numpy.array([0.0, 0.25, 0.80, 1.0])
                # z = numpy.array([ 
                #    [ 1.4 ,  6.5 ,  1.5 ,  1.8 ],
                #    [ 8.9 ,  7.3 ,  1.1 ,  1.09],
                #    [ 4.5 ,  9.2 ,  1.8 ,  1.2 ]])
                # # you have to set kx and ky small for this small example dataset
                # # 3 is more usual and is the default
                # # s=0 will ensure this interpolates.  s>0 will smooth the data
                # # you can also specify a bounding box outside the data limits
                # # if you want to extrapolate
                # sp = interpolate.RectBivariateSpline(x, y, z, kx=2, ky=2, s=0)
                #
                # sp([0.60], [0.25])  # array([[ 7.3]])
                # sp([0.25], [0.60])  # array([[ 2.66427408]])

           
                vl[0] = vector(fnn, fpn)
                vl[1] = vector(fnp, fpp)

                # NOTE: WE STORE ONES-BASED INDICES IN THE FACE_GRID_INDICES LIST
                #    TO MATCH CALCULIX NOTATION
                fgi = FACE_GRID_INDICES[face_num - 1]


                # GET PRESSURES FROM BLOCK FACES IN CASE THERE ARE PER-FACE PRESSURES

                # THE SCALE FACTORS STORED IN THE GRID'S norm_uvw MEMBER VARIABLE
                # ARE ALWAYS MEASURED FROM THE MINIMUM PARAMETER PLANE
                
                for e in ftel :
                    #print 'LOOKING AT ELEMENT =', e

                    face_pressure = 0.0
                    # AVERAGE PRESSURE FOR ALL GRIDS ON THE ELEMENT FACE
                    for i in fgi :
                        #print 'FACE GRID INDEX=',i
                        g = e.gl[i-1]
                        #print 'LOOKING AT ELEMENT FACE GRID -=>', g
                        #print 'G.NORM_UVW =', g.norm_uvw
                        
                        # WE DON'T ACTUALLY USE THE MESH GRID'S PRESSURE VALUE
                        # BUT SAMPLE THE PRESSURE OFF THE BLOCK FACE AT THE SAME
                        # PARAMETER VALUES AS WHERE THE GRID WAS DEFINED AT.
                        # THIS ALLOWS FOR DISCONTINUITY IN THE APPLIED SURFACE PRESSURE
                        # AT THE BLOCK EDGES.
                        # ALL THE OTHER PARAMETERS ARE AVERAGED THROUGH OUT THE VOLUME
                        # OF THE BLOCK
                        
                        scale_factor_1 = g.norm_uvw[ interp_dirs[0] ]
                        scale_factor_2 = g.norm_uvw[ interp_dirs[1] ]
                        #print 'SCALE FACTOR 1 =', scale_factor_1
                        #print 'SCALE FACTOR 2 =', scale_factor_2
                        pl0 = vl[0].point_from_scale( scale_factor_1 )
                        pl1 = vl[1].point_from_scale( scale_factor_1 )
                        vl[2] = vector(pl0, pl1)
                        pl = vl[2].point_from_scale( scale_factor_2 )
                        gp = pl.get_param('PRES', 0.0)
                        face_pressure = face_pressure + gp
                    pass

                    face_pressure = face_pressure / float( len(fgi) )

                    
                    ## #print 'ftel len = ', len(ftel)
                    ## #for e in ftel :
                    ## #    print e
                    ## #pass

                
                    ## # FROM CALCULIX MANUAL...
                    ## #  face 1: 1-2-3-4 -> FTAGS.MIN_W - 1
                    ## #  face 2: 5-8-7-6 -> FTAGS.MAX_W - 2
                    ## #  face 3: 1-5-6-2 -> FTAGS.MIN_V - 4
                    ## #  face 4: 2-6-7-3 -> FTAGS.MAX_U - 8
                    ## #  face 5: 3-7-8-4 -> FTAGS.MAX_V - 16
                    ## #  face 6: 4-8-5-1 -> FTAGS.MIN_U - 32

                    if( b.defined_by_tag == BTAGS.EXPLODED ) :

                        if( ft == FTAGS.MAX_W ) :
                            if( face_pressure < 0.0 ) : continue # WILL BE APPLIED TO THE -W FACE
                        pass

                        if( ft == FTAGS.MIN_W ) :
                            if( face_pressure > 0.0 ) : continue # WILL BE APPLIED TO THE +W FACE
                            if( face_pressure < 0.0 ) :
                                face_pressure = face_pressure * (-1.0)
                            pass
                        pass

                    pass
                    
                    if( abs(face_pressure) < constants.TOL ) : continue
                    s = str(e.id) + ', P' + str(face_num) + ', ' + str(face_pressure) + '\n'
                    fp.write(s)
                    #print s

                pass
            pass

        pass
       
        tcnt = 0
        #fp.write('*TEMPERATURE' + '\n')
        #fp.write('*BOUNDARY' + '\n')
        for b in self.blocks :
            m = b.mesh
            tgl = m.gl.values()
            tgl = filter( ( lambda g: g is not None ), tgl)
            tgl.sort(key=lambda g : g.id)        
            for g in tgl :
                if( g.id <= 0 ) : continue
                temp = g.get_param('TEMP', self.initial_temp)
                if( abs(temp - self.initial_temp) < constants.TOL ) : continue
                if( tcnt == 0 ) : fp.write('*TEMPERATURE' + '\n')
                tcnt += 1
                #if( temp is None ) : continue
                fp.write(str(g.id) + ', ' + str(temp) + '\n')
                #fp.write(str(g.id) + ', 11, 11, ' + str(temp) + '\n') 
            pass
        pass



        fp.write('*' * 80 + '\n')

        fp.write('*NODE PRINT, NSET=NALL' + '\n')
        fp.write('U' + '\n')
        fp.write('*NODE FILE, OUTPUT=3D' + '\n')
        fp.write('U' + '\n')
        fp.write('*EL PRINT, ELSET=EALL' + '\n')
        fp.write('S' + '\n')
        fp.write('*EL FILE, OUTPUT=3D' + '\n')
        fp.write('S' + '\n')
        fp.write('*END STEP' + '\n')
        fp.write('*' * 80 + '\n')


        fp.close()
    pass

    #----------------------------------------------------------------------------------------------------

    def write_nset(self, fp, name, gl) :
        fp.write('*NSET, NSET=' + name + '\n')
        nn = []
        for g in gl :
            nn.append(g.id)
        pass
        nn.sort()
        write_list_csv(fp, nn, 8) 
    pass

    #----------------------------------------------------------------------------------------------------

pass

##########################################################################################################
