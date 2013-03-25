import sys
import os
import re
import math
import string
#import types
import cPickle as pickle
#import shelve
#import errno
#import shutil
#import glob
import time
#import textwrap
#import random
import copy
from datetime import date
#from subprocess import *
#import csv
#import win32com.client
#import pythoncom
from xml.dom import minidom
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
#import scipy as sp
#from pylab import *
#import collections
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

## def get_grid_list_from_tags_gen(self, tag, igl = None) :
##     rgl = []
##     if( igl is not None ) :
##         rgl = igl
##         rgl = filter( ( lambda gg: gg is not None ), rgl)
##         rgl = filter( ( lambda gg: gg.tags & tag ), rgl )
##     pass
##     return(rgl)
## pass


##########################################################################################################

## def get_grid_list_from_not_tags_gen(self, tag, igl = None) :
##     rgl = []
##     if( igl is not None ) :
##         rgl = igl
##         rgl = filter( ( lambda gg: gg is not None ), rgl)
##         rgl = filter( ( lambda gg: not ( gg.tags & tag ) ), rgl )
##     pass
##     return(rgl)
## pass

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
# GRID IS NOW DERIVED FROM POINT

class grid(point) :

    def __init__ (self, gid = 0, x = 0.0, y = 0.0 , z = 0.0, rest = {}) :
        super(grid, self).__init__(x, y, z, rest)
        self.name = 'G' + str(gid)
        self.id = gid
        self.bc = set()
        self.iu = 0
        self.iv = 0
        self.iw = 0
        self.tags = 0
        self.partner_grids = set() # FOR JOINING MESHES
        self.mesh = None
        self.el = set()
        self.results = {}
        self.max_results = {}
        self.min_results = {}
        self.avg_results = {}
        #self.extra_grid = None
    pass

    #----------------------------------------------------------------------------------------------------

    def __del__ (self) :
        super(grid, self).__del__()
        del self.bc
        del self.el
        del self.partner_grids
        del self.results
        del self.max_results
        del self.min_results
        del self.avg_results
    pass

    #----------------------------------------------------------------------------------------------------

    def calc_avg_parameter_data_from_gl(self, gl) :
        #print '>>> TOP >>> ', here(), '  NAME = ',self.name
        
        super(grid, self).__init__(0.0, 0.0, 0.0)

        ng = len(gl)
        if( ng == 0 ) : return

        
        #print 'GRID:', self

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
#                     self.iu, self.iv, self.tags, ts)
        f = 'G: {0} ->  {1}, {2}, {3}    ( {4}, {5}, {6} ) BC:{7} ({8}) {9}'
        s = f.format(self.id, self.v[LOC.X], self.v[LOC.Y], self.v[LOC.Z],
                     self.iu, self.iv, self.iw, bc_stags, ts, mn)
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
    # THIS SHOULD REALLY ONLY RETURN A LIST OF ONE GRID - SEE model.stitch()
    #   - MAY RETURN A REFERENCE TO ITSELF
    def get_active_partner_grids(self) :
        pgrids = [ self ] + list(self.partner_grids)
        
        # agl = filter( ( lambda g : ( g.tags & GTAGS.MERGED )
        #                and ( not ( g.tags & GTAGS.REMOVED ) ) ), pgrids )
        
        agl = filter( ( lambda g : ( not ( g.tags & GTAGS.REMOVED ) ) ), pgrids )
        return(agl)
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
        #self.thick = 0.0
        #self.volume_fraction = 1.0
        #self.temp = 0.0
        #self.pres = 0.0
        self.mat_name = ""
        self.iu = 0
        self.iv = 0
        self.iw = 0
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
        
        self.lul = None
        self.lvl = None
        self.luh = None
        self.lvh = None

        self.model = None
        self.block = None
        #self.vol = None
        
        self.neu = 0
        self.nev = 0
        self.new = 0
        #self.eid_offset = 0
        #self.eid_max = -1
        #self.eid_min = -1
        self.ngu = 0
        self.ngv = 0
        self.ngw = 0
        #self.gid_offset = 0
        #self.gid_max = -1
        #self.gid_min = -1
        #self.gi_cl_u = 0
        #self.gi_cl_v = 0
        #self.el = []
        #self.gl = []
        self.el = {}
        self.gl = {}
        #self.parms = {}
        #self.tag = 0
        self.u_elsets = []
        self.v_elsets = []
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
        #del self.parms
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
        print 'NG =', b.ng
        print 'NE =', b.ne



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
                    g.tags = 0 # HANDLED BELOW FOR CLARITY
                    g.iu = iu
                    g.iv = iv
                    g.iw = iw
                    g.mesh = self
                    self.gl[iu, iv, iw] = g
                pass # W LOOP
            pass # V LOOP
        pass # U LOOP
    
        print 'GOING TO ASSIGN GRID TAGS'

        # GRID TAGS
        for iu in range(b.ng[DIR.U]) :
            
            for iv in range(b.ng[DIR.V]) :
            
                for iw in range(b.ng[DIR.W]) :

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


                    
                    self.gl[iu, iv, iw].tags = gtags
                    #print self.gl[iu, iv, iw]
                    
                pass # W LOOP
            pass # V LOOP
        pass # U LOOP

        # TBD: PROBABLY SHOULD WORK OUT NOT ALLOCATING THESE GRIDS
    
        #self.print_grids()
        #print 'DELETE UN-NEEDED GRIDS'
        dgl = self.get_grid_list_from_tags(GTAGS.PENDING_DELETE)
        for g in dgl :
            del self.gl[g.iu, g.iv, g.iw]
            self.gl[g.iu, g.iv, g.iw] = None
        pass
        del dgl
        #self.print_grids()
        
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

                    e.iu = iu
                    e.iv = iv
                    e.iw = iw
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

## // FYI: BIT TWIDDLING
## //   stat |= BIT; // TURNS IT ON IF ITS OFF - LEAVES ON IF IT IS ON - FORCES ON
## //   stat &= BIT; // LEAVES IT ON IF ITS ON - LEAVES OFF IF IT IS OFF - MASK TO SEE IF BIT IS ON 
## //   stat &= ~BIT; // TURNS IT OFF IF ITS ON - LEAVES OFF IF IT IS OFF - FORCES OFF
## //   stat ^= BIT; // TURNS IT ON IF ITS OFF - TURNS IT OFF ITS ON - TWIDDLE
## // NOT TO USEFUL ONES (HARD TO DESCRIBE TOO) :
## //   stat |= ~BIT; 
## //                      stat   BIT
## //                        0     0  -> 1
## //                        0     1  -> 0
## //                        1     0  -> 1
## //                        1     1  -> 1  
## //TURNS THEM ALL ON EXCEPT CURRENT BIT IF IT IS ON AND ORIGINAL WAS NOT
## //   stat ^= ~BIT; 
## //                      stat   BIT
## //                        0     0  -> 1
## //                        0     1  -> 0
## //                        1     0  -> 0
## //                        1     1  -> 1  

    def calc_grid_to_element_back_pointers(self) :
        for e in self.el.values() :
            for g in e.gl :
                g.el.add(e)
            pass
        pass
    pass

    #----------------------------------------------------------------------------------------------------

    def get_grid_from_indices(self, iu, iv, iw = 0, igl = None) :
        #print 'LOOKING FOR GRID AT INDEX =', iu, iv, iw

        rgl = self.gl;
        if( igl is not None ) :
            rgl = igl
        pass

        
        rgl = filter( ( lambda gg: gg.iu == iu and gg.iv == iv and gg.iw == iw ), rgl )

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
        #print 'LOOKING FOR GRID(S) AT INDEX =', iu, iv, iw
        #  REF: PYTHON library reference
        # Shallow copies of dictionaries can be made using dict.copy(),
        #   and of lists by assigning a slice of the entire list, for
        #   example, copied_list = original_list[:]        

        rgl = self.gl;
        if( igl is not None ) :
            rgl = igl
        pass
    
        if( iu is not None ) :
            rgl = filter( ( lambda gg: gg.iu == iu ), rgl )
        pass
        if( iv is not None ) :
            rgl = filter( ( lambda gg: gg.iv == iv ), rgl )
        pass
        if( iw is not None ) :
            rgl = filter( ( lambda gg: gg.iw == iw ), rgl )
        pass
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
            rel = filter( ( lambda ee: ee.iu == iu ), rel )
        pass
        if( iv is not None ) :
            rel = filter( ( lambda ee: ee.iv == iv ), rel )
        pass
        if( iw is not None ) :
            rel = filter( ( lambda ee: ee.iw == iw ), rel )
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
            #tgl = gg.get_active_grids()
            #if( len(tgl) > 0 ) :
            #    print ' '
            #    print gg
            #    print 'NACTIVE = ', len(tgl)
            #    print 'GRID ' + str(gg.id) + ' ACTIVE->', tgl
            #    print 'NPARTNERS = ', len(gg.partner_grids)
            #    print 'partner_grids = ', gg.partner_grids
            #pass
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
            #els = elset( ee.single_element_set_name, ee.id, self)
            #els.append(ee)
            #elsets_list.append( els )
            #fp.write('*ELEMENT, TYPE=S8R, ELSET=' + ee.single_element_set_name + '\n')
            fp.write('*ELEMENT, TYPE=C3D20R, ELSET=' + ee.single_element_set_name + '\n')
            #fp.write('*ELEMENT, TYPE=S8, ELSET=' + ee.single_element_set_name + '\n')
            ee.save_abq(fp)
        pass
        #return(elsets_list)
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
        self.face_tag = ref_block_tag
        self.hinge_dir = None

        # GRID LIST OF GRIDS ON HINGE LINE - FOUND NO MATTER WHAT KIND OF JOIN IT IS
        self.gl = []
        
        # GRID LIST OF GRIDS ON JOIN FACE - FOUND NO MATTER WHAT KIND OF JOIN IT IS
        self.face_gl = []

        # GRID LIST OF GRIDS ON LINE ON FACE ORTHOGONAL TO HINGE LINE
        self.ortho_gl = []

        self.edge_vector = np.zeros(3)
        self.ortho_vector = np.zeros(3)

        self.face_dirs = [-1, -1]

        self.get_face_directions()
        
        self.name = self.parent.name + '_' + self.block.name
    pass

    #----------------------------------------------------------------------------------------------------

    def get_face_directions(self) :
        self.face_dirs = [-1, -1]
        if( ( self.face_tag == FTAGS.MAX_U )
            or ( self.face_tag == FTAGS.MIN_U ) ) :
            self.face_dirs = [ DIR.V, DIR.W ]
        pass
        if( (  self.face_tag == FTAGS.MAX_V )
            or (  self.face_tag == FTAGS.MIN_V ) ) :
             self.face_dirs = [ DIR.U, DIR.W ]
        pass
        if( (  self.face_tag == FTAGS.MAX_W )
            or (  self.face_tag == FTAGS.MIN_W ) ) :
             self.face_dirs = [ DIR.U, DIR.V ]
        pass
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
        #self.block_a = ba
        #self.block_a_tag = ba_tag
        #self.block_b = bb
        #self.block_b_tag = bb_tag
        self.type = jtype

        self.master_tag = ''

        self.model = parent_model

        #self.jface_a.hinge_dir = None
        #self.jface_b.hinge_dir = None
        
        #self.a_gl = [] # GRID LIST OF "A" BLOCK GRIDS ON HINGE LINE - FOUND NO MATTER WHAT KIND OF JOIN IT IS
        #self.b_gl = [] # GRID LIST OF "B" BLOCK GRIDS ON HINGE LINE - FOUND NO MATTER WHAT KIND OF JOIN IT IS
        #self.jface_a.face_gl = [] # GRID LIST OF "A" BLOCK GRIDS ON JOIN FACE - FOUND NO MATTER WHAT KIND OF JOIN IT IS
        #self.jface_b.face_gl = [] # GRID LIST OF "B" BLOCK GRIDS ON JOIN FACE - FOUND NO MATTER WHAT KIND OF JOIN IT IS

        # THESE ARE SET EQUAL TO a_gl OR b_gl DEPENDING ON THIER LENGTH (master_gl IS THE LONGER)
        self.master = None 
        self.slave = None
        
        self.extra_gl = {} # THIS IS NOW A DICTIONARY WHOSE INDEX IS THE CORRESPONDING MASTER GRID OBJECT

        self.edge_vector = np.zeros(3)
        
    pass

    #----------------------------------------------------------------------------------------------------

    def __del__ (self) :
        del self.jface_a
        del self.jface_b
        #del self.jface_a.face_gl
        #del self.jface_b.face_gl
        del self.extra_gl
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
    def get_face_hinge_line_tags(block_face_tag, hinge_dir_tag) :

        hinge_grid_tag = -1
        ortho_grid_tag = -1
        ortho_grid_select_tag = -1

        # TBD - SHOULD SEE IF WE CAN "AND" THE ORTHO GRID TAGS  - SHOULD SPEED THINGS US 
        
        if( block_face_tag == GTAGS.MIN_U  ) :
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
    
        if( block_face_tag == GTAGS.MAX_U ) :
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
    
        if( block_face_tag == GTAGS.MIN_V ) :
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
    
        if( block_face_tag == GTAGS.MAX_V ) :
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
    
        if( block_face_tag == GTAGS.MIN_W ) :
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
    
        if( block_face_tag == GTAGS.MAX_W ) :
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

    # POPULATE THE GRID LISTS FOR JOIN FACES/LINES
    def get_grids(self) :

        a_ortho_grid_tag = -1
        b_ortho_grid_tag = -1
        a_ortho_grid_select_tag = -1
        b_ortho_grid_select_tag = -1

        
        ma = self.jface_a.block.mesh
        mb = self.jface_b.block.mesh

        print 'NUMBER OF GRIDS IN MESH A =', len(ma.gl)
        print 'NUMBER OF GRIDS IN MESH B =', len(mb.gl)

        print 'A FACE TAG = ', self.jface_a.face_tag
        print 'B FACE TAG = ', self.jface_b.face_tag


        self.jface_a.face_gl = ma.get_grid_list_from_tags(self.jface_a.face_tag)
        self.jface_b.face_gl = mb.get_grid_list_from_tags(self.jface_b.face_tag)



        
        if( self.type == JTAGS.HINGE ) :
        
            a_hinge_grid_tag, a_ortho_grid_tag, a_ortho_grid_select_tag = join.get_face_hinge_line_tags(self.jface_a.face_tag, self.jface_a.hinge_dir)
            self.jface_a.gl = ma.get_grid_list_from_tags(a_hinge_grid_tag, self.jface_a.face_gl)
            self.jface_a.ortho_gl =  ma.get_grid_list_from_tags(a_ortho_grid_tag, self.jface_a.face_gl)


            b_hinge_grid_tag, b_ortho_grid_tag, b_ortho_grid_select_tag = join.get_face_hinge_line_tags(self.jface_b.face_tag, self.jface_b.hinge_dir)
            self.jface_b.gl = mb.get_grid_list_from_tags(b_hinge_grid_tag, self.jface_b.face_gl)
            self.jface_b.ortho_gl =  ma.get_grid_list_from_tags(b_ortho_grid_tag, self.jface_b.face_gl)

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


            print 'NUMBER OF HINGE FACE GRIDS IN MESH A =', len(self.jface_a.face_gl)
            print 'NUMBER OF HINGE FACE GRIDS IN MESH B =', len(self.jface_b.face_gl)
            print 'NUMBER OF HINGE GRIDS IN MESH A =', len(self.jface_a.gl)
            print 'NUMBER OF HINGE GRIDS IN MESH B =', len(self.jface_b.gl)

        pass

        if( self.type == JTAGS.MERGE ) :
            
            for g in self.jface_a.face_gl : 
                g.tags |= GTAGS.MERGED_FACE
            pass

            for g in self.jface_b.face_gl : 
                g.tags |= GTAGS.MERGED_FACE
            pass
        
        pass
    
        ## print '_'*80
        ## print 'A GRIDS:', self.jface_a.face_tag, self.block_a.name
        ## c = -1
        ## for ag in self.jface_a.gl :
        ##    c += 1
        ##    print c, ag
        ## pass
    
        ## print
        
        ## print 'B GRIDS:', self.jface_b.face_tag, self.jface_b.block.name
        ## c = -1
        ## for bg in self.jface_b.gl :
        ##    c += 1
        ##    print c, bg
        ## pass

        ## print '_'*80
        return( a_ortho_grid_tag, b_ortho_grid_tag, a_ortho_grid_select_tag, b_ortho_grid_select_tag)
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
    # ASSOCIATES THE GRIDS ON THE MORE DENSE HINGE LINE (MASTER) WITH THE ONES ON THE LESS DENSE HINGE LINE (SLAVE)
    def find_partners(self) :
        print '>>> TOP >>> ', here()
        
        ma = self.jface_a.block.mesh
        mb = self.jface_b.block.mesh
        
        a_ortho_grid_tag, b_ortho_grid_tag, a_ortho_grid_select_tag, b_ortho_grid_select_tag = self.get_grids()

        print '_'*80
        print 'A FACE GRIDS:', self.jface_a.face_tag, self.jface_a.name
        c = -1
        for ag in self.jface_a.face_gl :
           c += 1
           print c, ag
        pass
    
        print
        
        print 'B FACE GRIDS:', self.jface_b.face_tag, self.jface_b.name
        c = -1
        for bg in self.jface_b.face_gl :
           c += 1
           print c, bg
        pass

        print '_'*80




        
        #    __ _______  ___________
        #   / // /  _/ |/ / ___/ __/
        #  / _  // //    / (_ / _/
        # /_//_/___/_/|_/\___/___/
        #

        if( self.type == JTAGS.HINGE ) :
            print 'HINGE JOIN'
            
            a_min_grid_tag = join.get_min_grid_tag(self.jface_a.hinge_dir)
            b_min_grid_tag = join.get_min_grid_tag(self.jface_b.hinge_dir)
            
            a_start_grid_list = ma.get_grid_list_from_tags(a_min_grid_tag, self.jface_a.gl)
            if( len(a_start_grid_list) != 1 ) :
                print 'ERROR:  COULD NOT FIND A UNIQUE END POINT FOR A HINGE ON BLOCK ', self.jface_a.block.name
                print '  Number of grids found = ', len(a_start_grid_list)
                return
            pass
            a_start_grid = a_start_grid_list[0]
            self.jface_a.gl.sort(key=lambda g : np.linalg.norm(g.v - a_start_grid.v))
            a_end_grid = self.jface_a.gl[-1]

            print 'A START GRID = ', a_start_grid
            print 'A END GRID = ', a_end_grid


            ## for g in self.jface_a.gl :
            ##     print 'LOOKING AT GRID:', g
            ##     dv = a_start_grid.v - g.v
            ##     print 'DV =', dv
            ##     print 'NORM=', np.linalg.norm(g.v - a_start_grid.v)
            ## pass
            
            self.jface_a.edge_vector = a_end_grid.v - a_start_grid.v
            
            a_ortho_grid_list = ma.get_grid_list_from_tags( a_ortho_grid_select_tag, self.jface_a.ortho_gl )
            if( len(a_ortho_grid_list) != 1 ) :
                print 'ERROR:  COULD NOT FIND A UNIQUE ORTHO GRID POINT DO DEFINE FACE COORDINATE SYSTEM ', self.jface_a.name
                print '  Number of grids found = ', len(a_ortho_grid_list)
                return
            pass                
            a_ortho_grid = a_ortho_grid_list[0]    

            print 'A SORTED GRIDS ='
            for g in self.jface_a.gl :
                print g
            pass

            print 'A ORTHO GRID = ', a_ortho_grid
            print 'A EDGE VECTOR = ', self.jface_a.edge_vector            

            self.jface_a.ortho_vector = a_ortho_grid.v - a_start_grid.v
            print 'A ORTHO VECTOR = ', self.jface_a.ortho_vector            

            ###--

            
            b_start_grid_list = ma.get_grid_list_from_tags(b_min_grid_tag, self.jface_b.gl)
            if( len(b_start_grid_list) != 1 ) :
                print 'ERROR:  COULD NOT FIND A UNIQUE END POINT FOR A HINGE ON BLOCK ', self.jface_b.block.name
                print '  Number of grids found = ', len(b_start_grid_list)
                return
            pass
            b_start_grid = b_start_grid_list[0]
            self.jface_b.gl.sort(key=lambda g : np.linalg.norm(b_start_grid.v - g.v))
            b_end_grid = self.jface_b.gl[-1]

            print 'B START GRID = ', b_start_grid
            print 'B END GRID = ', b_end_grid
            
            self.jface_b.edge_vector = b_end_grid.v - b_start_grid.v

            b_ortho_grid_list = mb.get_grid_list_from_tags( b_ortho_grid_select_tag, self.jface_b.ortho_gl )
            if( len(b_ortho_grid_list) != 1 ) :
                print 'ERROR:  COULD NOT FIND A UNIQUE ORTHO GRID POINT DO DEFINE FACE COORDINATE SYSTEM ', self.jface_b.name
                print '  Number of grids found = ', len(b_ortho_grid_list)
                return
            pass                
            b_ortho_grid = b_ortho_grid_list[0]    

            print 'B SORTED GRIDS ='
            for g in self.jface_b.gl :
                print g
            pass

            print 'B ORTHO GRID = ', b_ortho_grid
            print 'B EDGE VECTOR = ', self.jface_b.edge_vector            

            self.jface_b.ortho_vector = b_ortho_grid.v - b_start_grid.v
            print 'B ORTHO VECTOR = ', self.jface_b.ortho_vector            

            # WE REALLY HAVE TERMS MASTER AND SLAVE BACKWARDS
            # THE MASTER IS CALLED THE MASTER BECAUSE IT HAS A DENSER MESH
            #  IN THE HINGE DIRECTION AND THEREFORE IT DETERMINES THE NUMBER
            # OF EQUARTIONS TO BE WRITTEN.  BUT THE MASTER'S GRIDS (ADD EXTRA GRIDS)
            # ARE REALLY THE DEPENDENT VARIABLES WHOSE DISPLACEMENTS ARE DETERMINED
            # BY THE SLAVES GRID DISPLACEMENTS

            # SET THE MASTER FACE TO BE THE ONE WITH THE MOST GRID POINTS ON THE HINGE LINE 
            if( len(self.jface_a.gl) >= len(self.jface_b.gl) ) :
                self.master_tag = 'A'
                self.master = self.jface_a
                self.slave = self.jface_b
            else :
                self.master_tag = 'B'
                self.master = self.jface_b
                self.slave = self.jface_a
            pass

            for mg in self.master.gl :
                mg.tags |= GTAGS.JOIN_MASTER
            pass
            for sg in self.slave.gl :
                sg.tags |= GTAGS.JOIN_SLAVE
            pass

            coincident_tol = 0.01 # 1% of hinge line length
            self.edge_vector = self.master.edge_vector
            print 'EDGE VECTOR = ', self.edge_vector

            working_slave_gl = copy.copy(self.slave.gl)
            for mi, mg in enumerate(self.master.gl) :
                working_slave_gl.sort(key=lambda g : np.linalg.norm(g.v - mg.v))
                pg = working_slave_gl[0]
                mg.partner_grids.add(pg)
                dv = np.linalg.norm(pg.v - mg.v)
                #print 'MASTER GRID = ', mg.id
                #print 'PARTNER GRID =', pg.id
                #print 'COINCIDENT TOL = ', np.linalg.norm(self.edge_vector) * coincident_tol
                #print ' DV = ', dv
                if( dv > ( np.linalg.norm(self.edge_vector) * coincident_tol ) ) :
                     #print '  -- NOT COINCIDENT'
                     pg = working_slave_gl[1]
                     mg.partner_grids.add(pg)
                     #print '   SECOND PARTNER GRID =', pg.id
                pass
                # NOT SURE WE NEED THE EXTRA GRIDS BUT WE WILL KEEP THEM FOR NOW
                eg = grid(-1 * mg.id, mg.v[DIR.U], mg.v[DIR.V], mg.v[DIR.W])
                eg.tags |= GTAGS.HINGED_EXTRA
                self.extra_gl[mg] = eg
            pass

            print '_'*80
            print 'A HINGE LINE GRIDS:', self.jface_a.face_tag, self.jface_a.name
            c = -1
            for ag in self.jface_a.gl :
               c += 1
               print c, ag
            pass

            print

            print 'B HINGE LINE GRIDS:', self.jface_b.face_tag, self.jface_b.name
            c = -1
            for bg in self.jface_b.gl :
               c += 1
               print c, bg
            pass

            print '_'*80

        pass  # END JTAG.HINGE


    #   __  __________  _________
    #  /  |/  / __/ _ \/ ___/ __/
    # / /|_/ / _// , _/ (_ / _/
    #/_/  /_/___/_/|_|\___/___/
    #




        if( self.type == JTAGS.MERGE ) :
            print 'MERGE JOIN'
            
            # SET THE MASTER FACE TO BE THE ONE WITH THE MOST GRID POINTS ON THE FACE
            if( len(self.jface_a.face_gl) >= len(self.jface_b.face_gl) ) :
                self.master_tag = 'A'
                self.master = self.jface_a
                self.slave = self.jface_b
            else :
                self.master_tag = 'B'
                self.master = self.jface_b
                self.slave = self.jface_a
            pass

            for mg in self.master.gl :
                mg.tags |= GTAGS.JOIN_MASTER
            pass
            for sg in self.slave.gl :
                sg.tags |= GTAGS.JOIN_SLAVE
            pass

            working_slave_gl = copy.copy(self.slave.face_gl)
            # REMOVE ELEMENT_MID_EDGE NODES - USE ONLY CORNER NODES ON THE FACE
            #  SO WE CAN GET WRITE WEIGHTED EQUATIONS IN BOTH PARAMETRIC DIRECTIONS
            working_slave_gl = self.slave.block.mesh.get_grid_list_from_tags(GTAGS.ELEMENT_CORNER, working_slave_gl)

            print 'CORNER SLAVE GRIDS...'
            for si, sg in enumerate(working_slave_gl) :
                print si, sg
            pass


            
            for mi, mg in enumerate(self.master.face_gl) :
                print '\nLOOKING AT MASTER GRID...', mg
                working_slave_gl.sort(key=lambda g : np.linalg.norm(g.v - mg.v))
                pg = working_slave_gl[0] # THIS IS THE NEAREST SLAVE GRID

                mg.partner_grids.add(pg)
                vv = pg.v - mg.v
                dv = np.linalg.norm(vv)

                print 'FOUND CLOSE SLAVE GRID :', pg


                # WE ARE SITTING ON A POINT - WE HAVE AN 1 TO 1 RELATIONSHIP
                if( dv < constants.TOL ) :
                    print 'NEARLY IDENTTICAL SLAVE POINT - GOING TO MAP 1 TO 1  dv = ', dv
                    continue
                pass


                pg_uvw =  np.zeros( (4, 3), dtype=np.int32 )
                pg_uvw[0] = [ pg.iu, pg.iv, pg.iw ]


                # TBD -- WE NEED TO CHECK IF WE ARE RIGHT ON THE EDGE OF AN ELEMENT 

                # LOOK AT NEIGHBORING GRIDS IN THE FACE'S 2-D PARAMETER SPACE
                #  AND FIND THE CLOSET ONES TO THE POINT OF THE MASTER FACE WE ARE CONSIDERING
                i = 0
                for dd in self.slave.face_dirs :
                    i = i + 1
                    uvw = [None, None, None]
                    uvw[dd] = pg_uvw[0][dd]
                    line_o_grids = self.slave.block.mesh.get_grid_list_from_indices(
                        uvw[DIR.U], uvw[DIR.V], uvw[DIR.W], working_slave_gl)

                    if( pg in line_o_grids ) :
                        line_o_grids.remove(pg)
                    pass

                    line_o_grids.sort(key=lambda g : np.linalg.norm(g.v - pg.v))
                    npg = line_o_grids[0]
                    pg_uvw[i] =  [ npg.iu, npg.iv, npg.iw ]
                    mg.partner_grids.add(npg)
                pass

                ## print 'PARTNER GRIDS... (#=', len(mg.partner_grids)
                ## for tpg in mg.partner_grids :
                ##     print tpg
                ##     ell = list(tpg.el)
                ##     print ' BELONGS TO ELEMENTS:'
                ##     for ee in ell :
                ##         print '  ', ee.id,
                ##     pass
                ##     print
                ## pass


                ## tes = pg.el
                ## for tpg in mg.partner_grids :
                ##     tes = tes.intersection(tpg.el)
                ## pass
                ## print 'COMMON ELEMENTS =  (#=',len(tes)
                ## tel = list(tes)
                ## for ee in tel :
                ##     print '  ', ee.id,
                ## pass
                ## print
            

                # MAKE SURE THEY ALL BELONG TO THE SAME ELEMENT
                tes = pg.el
                for tpg in mg.partner_grids :
                    tes = tes.intersection(tpg.el)
                pass

                if( len(tes) == 1 ) :
                    if( len(mg.partner_grids) == 3 ) :
                        # FIND THE DIAGONAL GRID ON THE ELEMENT FACE 
                        duvw1 = pg_uvw[1] - pg_uvw[0]
                        duvw2 = pg_uvw[2] - pg_uvw[0]
                        #print 'DUVW1 = ', duvw1
                        #print 'DUVW2 = ', duvw2
                        pg_uvw[3] = pg_uvw[0] + duvw1 + duvw2
                        #print pg_uvw
                        npg = self.slave.block.mesh.gl[ pg_uvw[3][0], pg_uvw[3][1], pg_uvw[3][2] ]
                        mg.partner_grids.add(npg)
                        print 'DIAGONAL PARTNER...',
                        print npg

                        vv = mg.v - mg.partner_grids[0].v
                        vv1 =  mg.partner_grids[1].v - mg.partner_grids[0].v
                        vv1_uv = 
                        vv2 =  mg.partner_grids[3].v - mg.partner_grids[2].v



                        
                    pass
                pass

                print 'PARTNER GRIDS... (#=', len(mg.partner_grids)
                for tpg in mg.partner_grids :
                    print tpg
                    ell = list(tpg.el)
                    print ' BELONGS TO ELEMENTS:'
                    for ee in ell :
                        print '  ', ee.id,
                    pass
                    print
                pass

                tes = pg.el
                for tpg in mg.partner_grids :
                    tes = tes.intersection(tpg.el)
                pass
                print 'COMMON ELEMENTS =  (#=',len(tes)
                tel = list(tes)
                for ee in tel :
                    print '  ', ee.id,
                pass
                print
                
                



                ## print 'COMMON'
                ## comcome = set(tel)
                ## tel = list(comcome)
                ## for ee in tel :
                ##     print '  ', ee.id,
                ## pass
                ## print
               

                

                ## tgl = list(mg.partner_grids)
                ## com_el = tgl[0].el
                ## for tpg in mg.partner_grids :
                ##     com_el = com_el.intersection(tpg.el)
                ## pass
                ## print 'COMMON ELEMENTS = ', com_el

            pass


                ## deltas = [ 0, 0, 0 ]
                ## for i in self.slave.face_dirs :
                ##     deltas[i] = 1
                ## pass

                ## neighbor_indices = []
                ## if( ( self.face_tag == FTAGS.MAX_U )
                ##     or ( self.face_tag == FTAGS.MIN_U ) ) :
                ##     neighbor_indices.append( pg.iu, pg.iv + 1, pg.iw )  
                ##     neighbor_indices.append( pg.iu, pg.iv, pg.iw + 1 )  
                ##     neighbor_indices.append( pg.iu, pg.iv - 1, pg.iw )  
                ##     neighbor_indices.append( pg.iu, pg.iv, pg.iw - 1 )
                ## pass
                ## if( ( self.face_tag == FTAGS.MAX_V )
                ##     or ( self.face_tag == FTAGS.MIN_V ) ) :
                ##     neighbor_indices.append( pg.iu + 1, pg.iv, pg.iw )  
                ##     neighbor_indices.append( pg.iu, pg.iv, pg.iw + 1 )  
                ##     neighbor_indices.append( pg.iu - 1, pg.iv, pg.iw )  
                ##     neighbor_indices.append( pg.iu, pg.iv, pg.iw - 1 )
                ## pass
                ## if( ( self.face_tag == FTAGS.MAX_W )
                ##     or ( self.face_tag == FTAGS.MIN_W ) ) :
                ##     neighbor_indices.append( pg.iu + 1, pg.iv, pg.iw )  
                ##     neighbor_indices.append( pg.iu, pg.iv + 1, pg.iw )  
                ##     neighbor_indices.append( pg.iu - 1, pg.iv, pg.iw )  
                ##     neighbor_indices.append( pg.iu, pg.iv - 1, pg.iw )
                ## pass


                ## dirs = [ DIR.U, DIR.V, DIR.W ]
                ## for n in neighbor_indices :
                ##     for d in dirs :
                ##         if( neighbor_indices


                # NEED TO DETERMINE WHAT ELEMENT IF FALLS IN
                # IF IT DOES NOT FALL INSIDE AN ELEMENT THEN WE RIGHT EQUATIONS BASED ON THE LAPPED AREA ???
        
        pass # END JTAGS.MERGE


        # DEPRECATED - ASSUMES MATCHED GRID BETWEEN A AND B
        ## bgl = copy.copy(self.jface_b.gl)
        ## for ag in self.jface_a.gl :
        ##     #print
        ##     min_dv = constants.BIG_REAL
        ##     closest_grid_index = -1
        ##     c = closest_grid_index
        ##     #print 'LOOKING AT A GRID :', ag
        ##     for bg in bgl :
        ##         #print 'COMPARING TO B GRID :', bg
        ##         c += 1
        ##         dv = np.subtract(bg.v, ag.v)
        ##         #dvlen = np.linalg.norm(dv)
        ##         # http://stackoverflow.com/questions/9171158/how-do-you-get-the-magnitude-of-a-vector-in-numpy
        ##         # NEED FOR SPEED
        ##         #print 'DV=', dv
        ##         dvlen = np.sqrt( dv.dot( dv ) )
        ##         #print 'DVLEN = ', dvlen
        ##         if( dvlen < min_dv ) :
        ##             closest_grid_index = c
        ##             min_dv = dvlen
        ##             #print '*'
        ##         pass
        ##     pass
        ##     if( closest_grid_index >= 0 ) :
        ##         #print 'MIN DV = ', min_dv
        ##         #print 'A:', ag
        ##         #print 'B:', self.jface_b.gl[closest_grid_index]
         ##         # REMOVE FROM LOCAL LIST SO WE DON'T LOOK AT IT AGAIN
        ##         bg = bgl.pop(closest_grid_index)
        ##         if( self.type == JTAGS.MERGE ) :
        ##             ag.tags |= GTAGS.MERGED
        ##             bg.tags |= GTAGS.MERGED
        ##         pass
        ##         if( self.type == JTAGS.HINGE ) :
        ##             ag.tags |= GTAGS.HINGED
        ##             bg.tags |= GTAGS.HINGED
        ##         pass
        ##         #ag.partner_grids.append(bg)
        ##         #bg.partner_grids.append(ag)
        ##         ag.partner_grids.add(bg)
        ##         bg.partner_grids.add(ag)
        ##     pass
        ## pass
        
    pass

    #----------------------------------------------------------------------------------------------------

    def sew(self, join_type = None) :
        
        if( join_type != self.type ) :
            return
        pass
    
        if( self.type == JTAGS.MERGE ) :
            ggl = self.jface_a.gl + self.jface_b.gl

            ## for g in self.jface_a.gl :
            ##     print 'AGRID:', g
            ## pass
        
            ## print 

            ## for g in self.jface_b.gl :
            ##     print 'BGRID:', g
            ## pass

            print
        
            for g in ggl :
                pgrids = [ g ] + list(g.partner_grids)
                # FOR MERGED KEEP THE PARTNER (EQUIVALENT) GRID WITH THE LOWEST ID NUMBER
                # TAG THE REST AS REMOVED
                print
                print 'TOP: LOOKING AT GRID ', g.id
                print 'partner_grids = ', len( g.partner_grids )
                for gg in g.partner_grids :
                    print '    ', gg
                pass
                mg = min(pgrids, key = lambda pg : abs(pg.id) )
                print 'MIN GRID:', mg
                print 'MIN ID GRID partner_grids = ', len( mg.partner_grids )
                for gg in mg.partner_grids :
                    print '        ', gg
                pass
                agl = mg.get_active_partner_grids()
                print 'MIN ID GRID active_partner_grids = ', len( agl )
                for gg in agl :
                    print '   ', gg
                pass

            
                print 'LENGTH OF PGRIDS = ', len(pgrids)
                
                for pg in pgrids :
                    print '\n  >> LOOKING AT PGRID: ', pg
                    pg.tags |= GTAGS.MERGED
                    if( pg == mg ) :
                        print ' THIS IS THE MIN ID ONE - LETS  KEEP  IT\n'
                        continue
                    pass

                    print ' TAGGING AS GRID AS REMOVED: ', pg
                    pg.tags |= GTAGS.REMOVED
                    print '    NOW HAS BEEN TAGEGED AS REMOVED: ', pg

                    msh = pg.mesh
                    if( msh is None ) : continue

                    # MAKE SUBSTITUTIONS IN ELEMENT DEFINITIONS
                    pgel = msh.get_elements_using_grid(pg)
                    for ee in pgel :
                        for i in range(len(ee.gl)) :
                            tg = ee.gl[i]
                            if( tg == pg ) :
                                # SUB A REF TO THE GRID WE ARE KEEPING
                                print '\nSUB >>> ELEMENT ', ee.id, '  SUBBING GRID ', mg.id, 'FOR GRID ', ee.gl[i].id , ' (', tg.id, ')'
                                print '<<<ELEMENT WAS:', ee
                                ee.gl[i] = mg
                                print '>>>ELEMENT NOW', ee
                                pass
                            pass
                        pass
                 
                    pass
                pass

                print '\n***ALL ASSOCIATED GRIDS FOR FOR THE CURRENT GRID:', g
                for pg in pgrids :
                    c = ' '
                    if( not ( pg.tags & GTAGS.REMOVED ) ) :
                        c = '*'
                    pass
                    print c, '--->', pg
                pass
                print
            
            pass
        pass
    
    pass
pass

##########################################################################################################

# MODELS CAN HAVE MORE THAN ONE MESHES
class model(object) :

    def __init__ (self) :
        self.blocks = []
        self.meshes = []
        self.joins = []
        #  THESE ARE EXTRA GRIDS FOR HINGES
        #self.extra_gl = {} # THIS IS NOW A DICTIONARY WHOSE INDEX IS THE GRID WITH THE MINIMUM ID OF THE PARTNER_GRIDS
        self.materials = []
        self.initial_temp = constants.DEFAULT_INITIAL_TEMP;
        self.reference_temp = constants.DEFAULT_REFERENCE_TEMP; # SEE CALCULIX *EXPANSION CARD
        self.NPSS_materials = {};
        self.results_elements = [] # THIS IS A LIST OF ELEMENT CREATED WHEN THE RESULTS ARE READ IF THEY ARE NOT FOUND IN THE MESHES
        #self.stress_maxs = {}
        #self.stress_mins = {}
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

        maxs = [-1.0 * constants.BIG_REAL, -1.0 * constants.BIG_REAL, -1.0 * constants.BIG_REAL]
        mins = [ constants.BIG_REAL, constants.BIG_REAL, constants.BIG_REAL]
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
       

        print 'MAXS =', maxs
        print 'MINS =', mins
        
        ax.legend()
        #ax.set_xlim3d(mins[0], maxs[0])
        #ax.set_ylim3d(mins[1], maxs[1])
        #ax.set_zlim3d(mins[2], maxs[2])
        lo = min(mins) - 0.1 * abs(min(mins))
        hi = max(maxs) + 0.1 * abs(max(maxs))
        ax.set_xlim3d(lo, hi)
        ax.set_ylim3d(lo, hi)
        ax.set_zlim3d(lo, hi)

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

    ## def add_mesh(self, mname) :
    ##     m = mesh(mname)
    ##     m.model = self
    ##     self.meshes.append(m)
    ##     return(m)
    ## pass

    #----------------------------------------------------------------------------------------------------

    def add_material(self, mname) :
        m = material(mname)
        self.materials.append(m)
        return(m)
    pass

    #----------------------------------------------------------------------------------------------------

    ## def add_join(self, jname, ba, ba_tag, bb, bb_tag, jtype) :
    ##     j = join(jname, ba, ba_tag, bb, bb_tag, jtype, self)
    ##     self.joins.append(j)
    ##     return(j)
    ## pass

    #----------------------------------------------------------------------------------------------------

    def add_block_face_hinge_join(self, jname, ba, ba_face_tag, ba_hinge_dir, bb, bb_face_tag, bb_hinge_dir) :
        j = join(jname, ba, ba_face_tag, bb, bb_face_tag, JTAGS.HINGE, self)
        j.set_a_hinge_direction(ba_hinge_dir)
        j.set_b_hinge_direction(bb_hinge_dir)
        self.joins.append(j)
        return(j)
    pass

    #----------------------------------------------------------------------------------------------------


    def add_block_face_rigid_join(self, jname, ba, ba_face_tag, bb, bb_face_tag ) :
        j = join(jname, ba, ba_face_tag, bb, bb_face_tag, JTAGS.MERGE, self)
        j.set_a_hinge_direction(None)
        j.set_b_hinge_direction(None)
        self.joins.append(j)
        return(j)
    pass
        
    



    #----------------------------------------------------------------------------------------------------
    ## TBD
    ## def add_block_edge_hinge_join(self, ba, ba_edge_index, bb, bb_edge_index) :
    ##     j = join(ba, None, bb, None, JTAGS.HINGE, self)
    ##     j.set_a_hinge_direction(ba_edge_index)
    ##     j.set_b_hinge_direction(bb_edge_index)
    ##     self.joins.append(j)
    ##     return(j)
    ## pass

    #----------------------------------------------------------------------------------------------------

    ## def add_join(self, ma, ma_tag, mb, mb_tag, jtype) :
    ##     j = join(ma, ma_tag, mb, mb_tag, jtype, self)
    ##     self.joins.append(j)
    ##     return(j)
    ## pass

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
                #cards = dict.fromkeys(NPSS_required_fields, 0.0)
                #count_cards = {}
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

                # NPSS_required_fields = ['E', 'alpha', 'poissonRatio', 'rho']

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

    def read_stress_results_from_dot_frd_HEX(self, frd_fname) :

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


        # RIGHT NOW WE ARE ASSUMING WE ARE READING THE REUSLTS OF THE MESH THAT IS IN MEMORY
        #  IN THE FUTURE WE SHOULD BE ABLE TO READ THE FULL MESH IN FROM THE FRD FILE
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
            mosu = matp.props['SU'] / ( g.avg_results['VM'] * self.safety_factor_ultimate ) - 1.0
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
            mosu_mins = matp.props['SU'] / ( g.min_results['VM'] * self.safety_factor_ultimate ) - 1.0
            mosu_maxs = matp.props['SU'] / ( g.max_results['VM'] * self.safety_factor_ultimate ) - 1.0
            
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
            mosy = matp.props['SY'] / ( g.avg_results['VM'] * self.safety_factor_yield ) - 1.0
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
            mosy_mins = matp.props['SY'] / ( g.min_results['VM'] * self.safety_factor_yield ) - 1.0
            mosy_maxs = matp.props['SY'] / ( g.max_results['VM'] * self.safety_factor_yield ) - 1.0
  
        
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

    def stitch_HEX(self) :

        #max_gid = 0
        #for m in self.meshes :
        #    max_gid = max(max_gid, m.gid_max)
        #pass

                          
        for j in self.joins :
            print 'JOIN : ', j
            j.find_partners()
        pass

        # LIST OF ALL GRIDS
        #ggl = []
        #for m in self.meshes :
        #    ggl = ggl + m.gl
        #pass

        #print 'ELEMENTS BEFORE SEW...' 
        #for m in self.meshes :
        #    for e in m.el :
        #        print e
        #    pass
        #pass


        ## TBD - TBD - TBD
        ## # HAVE TO DO MERGES BEFORE WE DO HINGES 
        ## for j in self.joins :
        ##     j.sew(JTAGS.MERGE) # TBD
        ## pass
        ##
        ## # SWAP GRIDS TAGGED AS REMOVED WITH THE ACTIVE GRID AFTER ALL THE MERGE TYPE JOINS HAVE BEEN SEWN - THEN DO THE HINGES
        ## # AFTER THIS THE JOIN GRID LISTS (A & B) WILL ONLY CONTAIN GRIDS THAT ARE NOT TAGGED AS REMOVED
        ## # WE SHOULD NOT NEED TO CHECK FOR REMOVED GRIDS IN THE OUTPUT ROUTINES
        ## for  j in self.joins :
        ##     agllen = len(j.jface_a.gl)
        ##     for i in range(agllen) :
        ##         g = j.jface_a.gl[i]
        ##         if( g.tags & GTAGS.REMOVED ) :
        ##             tgl = g.get_active_partner_grids()
        ##             if( len(tgl) != 1 ) :
        ##                 print 'WARNING !!!  - A tgl is the wrong size -', len(tgl), ' for g=', g
        ##                 print here()
                        
        ##                 print 'PARTNER GRIDS = ', len(g.partner_grids) 
        ##                 cnt = 0
        ##                 for ggg in list(g.partner_grids) :
        ##                     print cnt, ggg
        ##                     cnt = cnt + 1
        ##                 pass
        ##                 print 
        ##             pass
        ##             #if( len( tgl ) == 0 ) : continue
        ##             g = tgl[0] # ASSUME THERE IS ONLY ONE - WHICH IS THE CASE IF THE REMOVE IS RESULTING FROM A MERGE
        ##             j.jface_a.gl[i] = g
        ##         pass
        ##     pass
        
        ##     bgllen = len(j.jface_b.gl)
        ##     for i in range(bgllen) :
        ##         g = j.jface_b.gl[i]
        ##         if( g.tags & GTAGS.REMOVED ) :
        ##             tgl = g.get_active_partner_grids()
        ##             if( len(tgl) != 1 ) :
        ##                 print 'WARNING !!!  - B tgl is the wrong size -', len(tgl), ' for g=', g
        ##                 print here()
                        
        ##                 print 'PARTNER GRIDS = ', len(g.partner_grids) 
        ##                 cnt = 0
        ##                 for ggg in list(g.partner_grids) :
        ##                     print cnt, ggg
        ##                     cnt = cnt + 1
        ##                 pass
                        
        ##             pass
        ##             #if( len( tgl ) == 0 ) : continue
        ##             g = tgl[0] # ASSUME THERE IS ONLY ONE - WHICH IS THE CASE IF THE REMOVE IS RESULTING FROM A MERGE
        ##             j.jface_b.gl[i] = g
        ##         pass
        ##     pass
        ## pass


        # HAVE TO PROCESS THE HINGES AFTER THE ALL THE MERGES  - IN CASE THERE IS SOME COMMONALITY IN THE GRIDS BEING AFFECTED
        for j in self.joins :
            j.sew(JTAGS.HINGE)
        pass

        gid = self.gid_max
        v = np.array([0.0, 0.0, 0.0])
        for j in self.joins :
            for kg, eg in j.extra_gl.items() :
                gid += 1
                print 'RE-ASSIGNING GRID ID ', gid, 'TO EXTRA GRID WHOSE KEY GRID IS ', kg.id
                eg.id = gid
            pass
        pass
        self.gid_max = gid

        #print 'ELEMENTS AFTER SEW...' 
        #for m in self.meshes :
        #    for e in m.el :
        #        print e
        #    pass
        #pass

        for b in self.blocks :
            b.mesh.calc_grid_to_element_back_pointers()
        pass

    pass

    #----------------------------------------------------------------------------------------------------

    def save_abq_HEX(self, fname) :
        print '>>> TOP >>> ', here()

        self.stitch_HEX()

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
        
        #ggl = ggl + list(self.extra_gl.values())

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


        # TBD
        # TAG UNUSED GRIDS WITH NEGATIVE IDS
        ## for g in ggl :
        ##     if( g.tags & GTAGS.REMOVED ) :
        ##         g.id = -1 * abs(g.id)
        ##     pass
        ## pass

        for m in self.materials :
            m.mat_props.sort(key=lambda prp: prp.temperature)
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

        ## egl = []
        ## for j in self.joins :
        ##     for eg in j.extra_gl.values() :
        ##         egl.append(eg)
        ##     pass
        ## pass
        ## egl.sort(key=lambda g: g.id)
                
        ## if( len(egl) > 0 ) :
        ##     fp.write('*** <BEGIN> EXTRA NODES FOR HINGE EQUATIONS' + '\n')
        ##     fp.write('*NODE, NSET=GHING' + '\n')
        ##     for eg in egl :
        ##         eg.save_abq(fp)
        ##     pass
        ##     fp.write('*** <END> EXTRA NODES FOR HINGE EQUATIONS' + '\n')
        ## pass

        for j in self.joins :
            egl = j.extra_gl.values()
            if( len(egl) > 0 ) :
                fp.write('*** <BEGIN> EXTRA NODES FOR HINGE EQUATIONS' + '\n')
                s = j.name + '_X'
                fp.write('*NODE, NSET=' + s + '\n')
                for eg in egl :
                    eg.save_abq(fp)
                pass
                fp.write('*** <END> EXTRA NODES FOR HINGE EQUATIONS' + '\n')
            pass
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



        for j in self.joins :
            for jf in j.jfaces :
                fp.write('** <BEGIN> ' + '~-~' * 15 + ' ' + jf.name  +' JOIN NODE SET\n')
                fp.write('*NSET, NSET=' + jf.name + '\n')
                nn = []
                for g in jf.face_gl :
                    nn.append(g.id)
                    pass
                nn.sort()
                write_list_csv(fp, nn, 8)
                fp.write('** <END> ' + '~-~' * 15 + ' ' + jf.name  +' JOIN NODE SET\n')
            pass
            fp.write('*NSET, NSET=' + j.name + '\n')
            s = ', '.join([ jf.name for jf in j.jfaces ])
            s = s + ', ' + j.name + '_X'
            fp.write(s + '\n')
        pass
    


        fp.write('*' * 80 + '\n')
        f = '{0}, {1}, {2}, {3}, {4}, {5}'

        for j in self.joins :
            jf = j.master
            fp.write('** <BEGIN> ' + '~-~' * 15 + ' ' + j.name  +' JOIN FACE COORDS\n')
            fp.write('*TRANSFORM, NSET=' + j.name + ', TYPE=R\n')
            s = f.format(jf.edge_vector[DIR.U], jf.edge_vector[DIR.V], jf.edge_vector[DIR.W],
                         jf.ortho_vector[DIR.U], jf.ortho_vector[DIR.V], jf.ortho_vector[DIR.W] )
            fp.write(s + '\n')
            fp.write('** <END> ' + '~-~' * 15 + ' ' + j.name  +' JOIN FACE COORDS\n')
        pass
    
        fp.write('*' * 80 + '\n')
        for j in self.joins :
            fp.write('** <BEGIN> ' + '~-~' * 5 + ' ' + jf.name  +' HINGE LINE EQUATIONS\n')
            for g in j.master.gl :
                pgl = list(g.partner_grids)
                lpg = len(g.partner_grids)
                print 'MASTER: GID =', g.id, ' NPARTNERS =', lpg
                eg = j.extra_gl[g]
                print 'EXTRA: GID =', eg.id
                
                if( lpg == 0 ) : continue
                
                if( lpg == 1 ) :
                    pg = pgl[0]
                    fp.write('*EQUATION\n')
                    for dof in [DOF.DX, DOF.DY, DOF.DZ] :
                        
                        fp.write( str(lpg+1) + '\n')
                        fp.write( str(eg.id) + ', ' + str(dof) + ', -1.0, '
                                  + str(pg.id) + ', ' + str(dof) + ', 1.0\n' )
                        
                        fp.write( '2\n')                        
                        fp.write( str(g.id) + ', ' + str(dof) + ', -1.0, '
                                  + str(eg.id) + ', ' + str(dof) + ', 1.0\n' )   
                    pass
                    ## for dof in [DOF.DX, DOF.DY, DOF.DZ] :
                    ##     fp.write( str(lpg+1) + '\n')
                    ##     fp.write( str(g.id) + ', ' + str(dof) + ', -1.0, '
                    ##               + str(pg.id) + ', ' + str(dof) + ', 1.0\n' )   
                    ## pass
                if( lpg == 2 ) :
                    pg1 = pgl[0]
                    pg2 = pgl[1]
                    fp.write('*EQUATION\n')
                    dv = np.linalg.norm(pg2.v - pg1.v)
                    dvg1 = np.linalg.norm(pg1.v - g.v)
                    dvg2 = np.linalg.norm(pg2.v - g.v)
                    f1 = 1.0 - dvg1 / dv
                    f2 = 1.0 - dvg2 / dv
                    for dof in [DOF.DX, DOF.DY, DOF.DZ] :

                        fp.write( str(lpg+1) + '\n')
                        fp.write( str(eg.id) + ', ' + str(dof)  + ', -1.0, '
                                  + str(pg1.id) + ', ' + str(dof) + ', ' + str(f1) + ', '
                                  + str(pg2.id) + ', ' + str(dof) + ', ' + str(f2) + ' \n' )
                        fp.write( '2\n')                        

                        fp.write( str(g.id) + ', ' + str(dof) + ', -1.0, '
                                  + str(eg.id) + ', ' + str(dof) + ', 1.0\n' )   
                    pass
                    ## for dof in [DOF.DX, DOF.DY, DOF.DZ] :
                    ##     fp.write( str(lpg+1) + '\n')
                    ##     fp.write( str(g.id) + ', ' + str(dof)  + ', -1.0, '
                    ##               + str(pg1.id) + ', ' + str(dof) + ', ' + str(f1) + ', '
                    ##               + str(pg2.id) + ', ' + str(dof) + ', ' + str(f2) + ' \n' )
                    ## pass
                    
                pass
            
                for pg in g.partner_grids :
                    print '    PARTNER GRID =', pg.id
                pass
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

        
        ## fp.write('*' * 80 + '\n')
        ## # TBD
        ## # WE MAY NEED A LOCAL COORDINATE SYSTEM HERE - TEST THIS !!!!
        ## f = '{agid}, {adof}, -1.0, {bgid}, {bdof}, 1.0'
        ## for j in self.joins :
        ##     for eg in j.extra_gl.values() :
        ##         fp.write('*EQUATION' + '\n')
        ##         for egp in eg.partner_grids :
        ##             # THIS ASSUMES HINGE TRANSFERS ONLY TRANSLATIONS IN X, Y AND Z - REALLY BALL AND SOCKET
        ##             for dof in [DOF.DX, DOF.DY, DOF.DZ] :
        ##                 # DON'T WRITE EQUATION REFERENCING FIXED GRID DOF
        ##                 if( ( dof in egp.bc ) or ( dof in eg.bc ) ) : continue
        ##                 s = f.format(agid = str(egp.id),
        ##                          adof = str(dof),
        ##                          bgid = str(eg.id),
        ##                          bdof = str(dof))
        ##                 fp.write('2' + '\n')
        ##                 fp.write(s + '\n')
        ##             pass           
        ##         pass
        ##     pass
        ## pass 

        if( self.do_nlgeom ) :
            fp.write('*STEP, INC=100, NLGEOM' + '\n')
            fp.write('*STATIC' + '\n')
            fp.write('0.1, 1.0' + '\n')
        else :
            fp.write('*STEP' + '\n')
            fp.write('*STATIC' + '\n')
        pass
    


        
        fp.write('*DLOAD' + '\n')

        ## for m in self.meshes :
        ##     for es in (m.u_elsets) :
        ##         pres = es.avg_grid.get_param('PRES')
        ##         if( pres is None ) : continue
        ##         fp.write(es.name + ', P, ' + str(pres) + '\n')
        ##     pass
        ## pass


        face_tags = [ FTAGS.MAX_U, FTAGS.MIN_U, FTAGS.MAX_V, FTAGS.MIN_V, FTAGS.MAX_W, FTAGS.MIN_W ]
        face_tags_exploded = [ FTAGS.MAX_W, FTAGS.MIN_W ]
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
                #print 'ft = ', ft
                face_num = int( ( math.log(ft) / math.log(2.0) ) + 1 )
                #print 'FACENUM = ', face_num

                ftel = filter( ( lambda ee : ee.tags & ft ), tel )
                
                #print 'ftel len = ', len(ftel)
                #for e in ftel :
                #    print e
                #pass

                
                # FROM CALCULIX MANUAL...
                #  face 1: 1-2-3-4 -> FTAGS.MIN_W - 1
                #  face 2: 5-8-7-6 -> FTAGS.MAX_W - 2
                #  face 3: 1-5-6-2 -> FTAGS.MIN_V - 4
                #  face 4: 2-6-7-3 -> FTAGS.MAX_U - 8
                #  face 5: 3-7-8-4 -> FTAGS.MAX_V - 16
                #  face 6: 4-8-5-1 -> FTAGS.MIN_U - 32

                
                fgi = FACE_GRID_INDICES[face_num - 1]

                for e in ftel :

                    facep = 0.0
                    for i in fgi :
                        # NOTE: WE STORE ONES-BASED INDICES IN THE FACE_GRID_INDICES LIST
                        #    TO MATCH CALCULIX NOTATION
                        g = e.gl[i-1] 
                        gp = g.get_param('PRES', 0.0)
                        #if( gp is None ) : continue
                        facep = facep + gp
                    pass

                    facep = facep / float( len(fgi) )

                    if( b.defined_by_tag == BTAGS.EXPLODED ) :

                        if( ft == FTAGS.MAX_W ) :
                            if( facep < 0.0 ) : continue # WILL BE APPLIED TO THE -W FACE
                        pass

                        if( ft == FTAGS.MIN_W ) :
                            if( facep > 0.0 ) : continue # WILL BE APPLIED TO THE +W FACE
                            if( facep < 0.0 ) :
                                facep = facep * (-1.0)
                            pass
                        pass

                    pass
                    if( abs(facep) < constants.TOL ) : continue
                    fp.write(str(e.id) + ', P' + str(face_num) + ', ' + str(facep) + '\n')

                pass
            pass
            

                
            ## for e in tel  :
            ##     pres = e.avg_grid.get_param('PRES')
            ##     if( pres is None ) : continue
            ##     if( abs(pres) < constants.TOL ) : continue
            ##     fp.write(str(e.id) + ', P, ' + str(pres) + '\n')
            ## pass
        pass
    
        ## for b in self.blocks :
        ##     m = b.mesh
        ##     tel = m.el.values()
        ##     tel.sort(key=lambda e : e.id)
        ##     for e in tel  :
        ##         pres = e.avg_grid.get_param('PRES')
        ##         if( pres is None ) : continue
        ##         if( abs(pres) < constants.TOL ) : continue
        ##         fp.write(str(e.id) + ', P, ' + str(pres) + '\n')
        ##     pass
        ## pass


        
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


        #  CALCULIX HAS A BUG - IT DOES NOT CALCULATE SHELL STRESSES RIGHT
        #  UNLESS YOU ASK FOR 'OUTPUT=3D'
        #    THIS WILL ONLY PUT THE CORRECT STRESSES IN THE FRD FILE FOR
        #    THE DYNAMICALLY CREATED 3D MESH
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

pass

##########################################################################################################
