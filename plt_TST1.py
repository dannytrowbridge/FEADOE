import doe.analysis as doe_analysis
import doe.gfea_model as fea
from doe.ggen_sup import point, vector
#import doe.constants
from doe.constants import LOC, DIR, DOF, GTAGS, ETAGS, FTAGS, JTAGS

import math

## http://code.activestate.com/recipes/52192-add-a-method-to-a-class-instance-at-runtime/

anal = doe_analysis.analysis('PLT_VT1')

#anal.add_indep_var('NEX', [5, 10, 15, 20, 25, 30])
anal.add_indep_var('NEX', [40])
#anal.add_indep_var('NEX', [1, 3, 5, 10, 20, 30, 40, 60, 80, 100 ])
#anal.add_indep_var('NEX', [25])
anal.add_indep_var('Pressure', [10000.0])
anal.add_indep_var('matty', ['Alloy_713LC'])

#####################################################################
                   
def define_model(self) :

    #self.do_nlgeom = True
    self.do_nlgeom = False
    self.set_initial_temp(0.0)

    #mat = self.add_material('MAT')
    #mat.set_number_of_output_divisions(6)
    #mat.add_props_at_temp(0.0, 1.44E-4, 35E6, 0.3, 5.9E-6)
    #mat.add_props_at_temp(50.0, 1.42E-4, 32E6, 0.3, 6.1E-6)
    #mat.add_props_at_temp(70.0, 1.39E-4, 30E6, 0.3, 6.16E-6)
    #mat.add_props_at_temp(100.0, 1.37E-4, 29E6, 0.3, 6.21E-6)

    mat = self.add_NPSS_material(self.matty)
    mat.set_number_of_output_divisions(1)

    #dmat = self.add_NPSS_material(self.matty)
    #dmat.set_number_of_output_divisions(7)

    # bottom left coords
    blx = 0.0
    bly = -1.0
    blz = 0.0

    dx = 20.0
    dy = 2.0
    dz = 0.0

    pp = {}

    th = 1.0

##     pp['000'] = point(blx + dx * 0.0, bly + dy * 0.0, blz + dz * 0.0,
##                       {'THICK': th, 'TEMP': 0.0, 'PRES': self.Pressure } )
##     pp['100'] = point(blx + dx * 1.0, bly + dy * 0.0, blz + dz * 0.0,
##                       {'THICK': th, 'TEMP': 0.0, 'PRES': self.Pressure } )
##     pp['010'] = point(blx + dx * 0.0, bly + dy * 1.0, blz + dz * 0.0,
##                       {'THICK': th, 'TEMP': 0.0, 'PRES': self.Pressure } )
## #    pp['001'] = point(blx + dx * 0.0, bly + dy * 0.0, blz + dz * 1.0,
## #                      {'THICK': 0.5, 'TEMP': 0.0, 'PRES': self.Pressure } )
##     pp['110'] = point(blx + dx * 1.0, bly + dy * 1.0, blz + dz * 0.0,
##                       {'THICK': th, 'TEMP': 0.0, 'PRES': self.Pressure } )
## #    pp['011'] = point(blx + dx * 0.0, bly + dy * 1.0, blz + dz * 1.0,
## #                      {'THICK': 0.5, 'TEMP': 0.0, 'PRES': self.Pressure } )
## #    pp['101'] = point(blx + dx * 1.0, bly + dy * 0.0, blz + dz * 1.0,
## #                      {'THICK': 0.5, 'TEMP': 0.0, 'PRES': self.Pressure } )
## #    pp['111'] = point(blx + dx * 1.0, bly + dy * 1.0, blz + dz * 1.0,
## #                      {'THICK': 0.5, 'TEMP': 0.0, 'PRES': self.Pressure } )

    common_params =  {'THICK': th, 'TEMP': 0.0, 'PRES': self.Pressure }



    ## ang = 30.0

    ## ca = math.cos(math.radians(ang))
    ## sa = math.sin(math.radians(ang))

    ## x = 0.0
    ## y = -1.0
    ## pp['000'] = point(x*ca - y*sa, x*sa + y*ca, 0.0, common_params)
    ## x = 20.0
    ## y = -1.0
    ## pp['100'] = point(x*ca - y*sa, x*sa + y*ca, 0.0, common_params)
    ## x = 0.0
    ## y = 1.0
    ## pp['010'] = point(x*ca - y*sa, x*sa + y*ca, 0.0, common_params)
    ## x = 20.0
    ## y = 1.0
    ## pp['110'] = point(x*ca - y*sa, x*sa + y*ca, 0.0, common_params)


    pp['000'] = point(0.0, 0.0, 0.0, common_params)
    pp['100'] = point(4.0, -2.0, 0.0, common_params)
    pp['010'] = point(0.0, 2.0, 0.0, common_params)
    pp['110'] = point(4.0, 5.0, 0.0, common_params)
    



    pp['001'] = point(0.0, 0.0, 3.0, {'THICK': th, 'TEMP': 0.0, 'PRES': self.Pressure } )
    pp['011'] = point(0.0, 3.0, 3.0, {'THICK': th, 'TEMP': 0.0, 'PRES': self.Pressure } )


    nex = self.NEX
    ney = 2
    eid_off = 0
    gid_off = 0

    ms = {}


    # EXTERNAL PRESSURE - NORMALS POINT IN
    
    mname = 'A'
    ms[mname] = self.add_block(mname)
    ms[mname].set_material(mat)
    ms[mname].define_block_from_mid_plane_points(pp['000'], pp['100'], pp['110'], pp['010'] )


    mname = 'B'
    ms[mname] = self.add_block(mname)
    ms[mname].set_material(mat)
    ms[mname].define_block_from_mid_plane_points(pp['000'], pp['010'], pp['011'], pp['001'] )

    self.plot_blocks()


    # BLOCK :  +/- 1 parametric CUBE CENTERED AT ZERO
    # EDGE INDICIES (su, sv, sw, eu, ev, ew) ALWAYS POINT FROM LO TO HI VALUE
    # FIRST THREE INDICIES ARE THE START POINT THE LAST THREE ARE THE END POINT INDICIES
    #  i.e.   EDGE(-1, 1, 1, 1, 1, 1) POINTS FROM POINT (-1, 1, 1) TO (1, 1, 1)
    #     THE EDGE VECTOR POINTS DOWN THE U AXIS ALONG THE EDGE WHERE V, ARE W ARE THIER MAXIMUM
    #  ONLY ONE OF THE INDICIES WILL CHANGE FROM (-) TO (+) FROM THE START POINT TO THE END POINT
    #  EACH BLOCK STORES 8 CORNER POINTS AND 12 EDGE VECTORS
    

    #                              BLOCK,   FACE_TAG,    HINGE_AXIS_DIR
    #self.add_block_face_join('AMB', ms['A'], FTAGS.MIN_U,
    #                               ms['B'], FTAGS.MIN_U,) 

    #                              BLOCK,   FACE_TAG,    HINGE_AXIS_DIR
    self.add_block_face_hinge_join('AHB', ms['A'], GTAGS.MIN_U, DIR.V,
                                   ms['B'], GTAGS.MIN_V, DIR.U) 

    # TBD
    # self.add_block_edge_hinge_join(ms['A'], (1, -1, -1, 1, 1, -1),
    #                                ms['B'], (-1, -1, -1, -1, 1, -1))

    self.generate_mesh()

    gl = ms['A'].get_grid_list_from_tags(GTAGS.MAX_U)
    gl = ms['A'].get_grid_list_from_tags(GTAGS.CENTER_W, gl)
    gl = gl + ms['B'].get_grid_list_from_tags(GTAGS.MAX_V)
    for g in gl :
         g.bc.add(DOF.DX)
         g.bc.add(DOF.DY)
         g.bc.add(DOF.DZ)
         #print 'CONSTRAINED GRID :', g
    pass
    #self.plot_mesh()
    
    print ' *** M O D E L   B U I L D   D O N E *** '

pass

#####################################################################

anal.analyze(define_model)



