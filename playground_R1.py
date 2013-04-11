import doe.analysis as doe_analysis
import doe.gfea_model as fea
from doe.ggen_sup import point, vector, face
#import doe.constants
from doe.constants import LOC, DIR, DOF, GTAGS, ETAGS, FTAGS, JTAGS

import math
import copy

## http://code.activestate.com/recipes/52192-add-a-method-to-a-class-instance-at-runtime/

anal = doe_analysis.analysis('PLAYGROUND')

#anal.add_indep_var('NEX', [5, 10, 15, 20, 25, 30])
anal.add_indep_var('NEX', [40])
#anal.add_indep_var('NEX', [1, 3, 5, 10, 20, 30, 40, 60, 80, 100 ])
#anal.add_indep_var('NEX', [25])
anal.add_indep_var('Pressure', [-10000.0])
anal.add_indep_var('matty', ['Alloy_713LC'])

#####################################################################
                   
def define_model(self) :

    #self.do_nlgeom = True
    self.do_nlgeom = False
    self.set_initial_temp(0.0)

    mat = self.add_NPSS_material(self.matty)
    mat.set_number_of_output_divisions(1)


    common_params =  {'TEMP': 0.0, 'PRES': 0.0 }
    common_params_L =  {'TEMP': 0.0, 'PRES': self.Pressure }

    # FROM CALCULIX MANUAL...
    #  face 1: 1-2-3-4 -> FTAGS.MIN_W - 1
    #  face 2: 5-8-7-6 -> FTAGS.MAX_W - 2
    #  face 3: 1-5-6-2 -> FTAGS.MIN_V - 4
    #  face 4: 2-6-7-3 -> FTAGS.MAX_U - 8
    #  face 5: 3-7-8-4 -> FTAGS.MAX_V - 16
    #  face 6: 4-8-5-1 -> FTAGS.MIN_U - 32

    # Must use faces (not points) so we can change pressures, etc. independently on each face  

##                            
##                                                   
##                8                        7    
##               .:........................  
##             .' |      F2             .'
##         5 .:....................._.-' |
##           |    |                 |6   |
##           |    |      F5(back)   |    |
##           |    |                 |    |                   
##           | F6 |    F3(front)    | F4 |
##           |    |                 |    |
##           |    |4                |    |
##           |   .'.................|....' 3
##           |  /        F1         | .'
##           `.:....................|'
##           1                      2
##                                  
## FACE NORMALS POINT IN
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

    p1 = point(0.0, 0.0, 0.0, common_params)
    p2 = point(5.0, 0.0, 0.0, common_params)
    p3 = point(5.0, 1.0, 0.0, common_params)
    p4 = point(0.0, 1.0, 0.0, common_params)
    p5 = point(0.0, 0.0, 1.0, common_params)
    p6 = point(5.0, 0.0, 1.0, common_params)
    p7 = point(5.0, 1.0, 1.0, common_params)
    p8 = point(0.0, 1.0, 1.0, common_params)



    p2_L = copy.deepcopy(p2)
    p2_L.set_param('PRES', self.Pressure)
    
    p3_L = copy.deepcopy(p3)
    p3_L.set_param('PRES', self.Pressure)

    p6_L = copy.deepcopy(p6)
    p6_L.set_param('PRES', self.Pressure)
    
    p7_L = copy.deepcopy(p7)
    p7_L.set_param('PRES', self.Pressure)



    
    # FROM CALCULIX MANUAL...
    #  face 1: 1-2-3-4 -> FTAGS.MIN_W - 1
    #  face 2: 5-8-7-6 -> FTAGS.MAX_W - 2
    #  face 3: 1-5-6-2 -> FTAGS.MIN_V - 4
    #  face 4: 2-6-7-3 -> FTAGS.MAX_U - 8
    #  face 5: 3-7-8-4 -> FTAGS.MAX_V - 16
    #  face 6: 4-8-5-1 -> FTAGS.MIN_U - 32

    f1 = face('F1')
    f1.define_face_from_points(p1, p2, p3, p4)
        
    f2 = face('F2')
    f2.define_face_from_points(p5, p8, p7, p6)

    f3 = face('F3')
    f3.define_face_from_points(p1, p5, p6, p2)
    
    #f4 = face('F4')
    #f4.define_face_from_points(p2, p6, p7, p3)
    
    f4_L = face('F4_L')
    f4_L.define_face_from_points(p2_L, p6_L, p7_L, p3_L)

    f5 = face('F5')
    f5.define_face_from_points(p3, p7, p8, p4)

    f6 = face('F6')
    f6.define_face_from_points(p4, p8, p5, p1)
    
    bname = 'A'
    blk_a = self.add_block(bname)
    blk_a.set_material(mat)
    blk_a.define_block_from_faces(f1, f2, f3, f4_L, f5, f6)

    self.plot_blocks()


    # BLOCK :  +/- 1 parametric CUBE CENTERED AT ZERO
    # EDGE INDICIES (su, sv, sw, eu, ev, ew) ALWAYS POINT FROM LO TO HI VALUE
    # FIRST THREE INDICIES ARE THE START POINT THE LAST THREE ARE THE END POINT INDICIES
    #  i.e.   EDGE(-1, 1, 1, 1, 1, 1) POINTS FROM POINT (-1, 1, 1) TO (1, 1, 1)
    #     THE EDGE VECTOR POINTS DOWN THE U AXIS ALONG THE EDGE WHERE V, ARE W ARE THIER MAXIMUM
    #  ONLY ONE OF THE INDICIES WILL CHANGE FROM (-) TO (+) FROM THE START POINT TO THE END POINT
    #  EACH BLOCK STORES 8 CORNER POINTS AND 12 EDGE VECTORS
    

    #                              NAME,  BLOCK,   GRID_TAG
    #self.add_block_face_rigid_join('ARB', ms['A'], GTAGS.MIN_U,
    #                                      ms['B'], GTAGS.MIN_V) 

    #                              BLOCK,   GRID_TAG,    HINGE_AXIS_DIR
    #self.add_block_face_hinge_join('AHB', ms['A'], GTAGS.MIN_U, DIR.V,
    #                               ms['B'], GTAGS.MIN_V, DIR.U) 

    # MUST GENERATE THE MESH BEFORE YOU CAN APPLY THE BOUNDARY CONDITIONS
    self.generate_mesh()

    #gla = ms['A'].get_grid_list_from_tags(GTAGS.MAX_U)
    #gla = ms['A'].get_grid_list_from_tags(GTAGS.CENTER_W, gla)
    
    glb = blk_a.get_grid_list_from_tags(GTAGS.MIN_U)
    glbw = blk_a.get_grid_list_from_tags(GTAGS.CENTER_W, glb)
    glbv = blk_a.get_grid_list_from_tags(GTAGS.CENTER_V, glb)

    print 'NUMBER OF GRIDS TO BE CONSTRAINED =', len(glb)



    for g in glb :
         g.bc.add(DOF.DX)
    pass


    for g in glbw :
         g.bc.add(DOF.DY)
    pass

    for g in glbv :
         g.bc.add(DOF.DZ)
    pass


    
    print ' *** M O D E L   B U I L D   D O N E *** '

pass

#####################################################################

anal.analyze(define_model)



