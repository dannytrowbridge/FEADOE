import sys
import doe.analysis as doe_analysis
import doe.gfea_model as fea
from doe.ggen_util import print_timing, enum, here, now, interpolate
from doe.ggen_sup import point, vector, face
from doe.constants import LOC, DIR, DOF, GTAGS, ETAGS, FTAGS, JTAGS

import math
import copy

anal = doe_analysis.analysis('NOZ_T2')
anal.add_indep_var('Pressure', [10000.0])
anal.add_indep_var('matty', ['Alloy_713LC'])

#####################################################################
@print_timing                   
def define_model(self) :

    #self.do_nlgeom = True
    self.do_nlgeom = False
    self.set_initial_temp(0.0)

    mat = self.add_NPSS_material(self.matty)
    mat.set_number_of_output_divisions(1)


    IN = 0
    OUT = 1

    common_params =  { 'PRES' : 0.0, 'VOL_FRAC' : 1.0, 'TEMP' : 0.0 }

    # FROM CALCULIX MANUAL...
    #  face 1: 1-2-3-4 -> FTAGS.MIN_W - 1
    #  face 2: 5-8-7-6 -> FTAGS.MAX_W - 2
    #  face 3: 1-5-6-2 -> FTAGS.MIN_V - 4
    #  face 4: 2-6-7-3 -> FTAGS.MAX_U - 8
    #  face 5: 3-7-8-4 -> FTAGS.MAX_V - 16
    #  face 6: 4-8-5-1 -> FTAGS.MIN_U - 32

    # Must use faces (not points) so we can change pressures, etc. independently on each face  

    txs = 0.0
    txe = 21.0

    tys = -2.0
    tye = 2.0

    tzs = 10.0
    tze = 11.0
    
    thref = 1.0

    hinge_x = 9.0
    hinge_z = 5.0

    ring_x = 17.0
    ring_z = interpolate(hinge_z, hinge_x, tze, txe, ring_x)

       
#  ---- ---- ---- ---- A BLOCK

    sx = [ txs, txs ]
    ex = [ hinge_x, hinge_x ]

    sy = [ tys, tys ]
    ey = [ tye, tye ]

    # ALWAYS ADD THE Z THICKNESS TO THE IN SURFACE TO GET THE OUT SURFACE
    sz = [ tzs, tzs + thref ]
    ez = [ hinge_z, hinge_z + thref]

    
    p1 = point(sx[IN], sy[IN], sz[IN], common_params)
    p2 = point(ex[IN], sy[IN], ez[IN], common_params)
    p3 = point(ex[IN], ey[IN], ez[IN], common_params)
    p4 = point(sx[IN], ey[IN], sz[IN], common_params)
    p5 = point(sx[OUT], sy[OUT], sz[OUT], common_params)
    p6 = point(ex[OUT], sy[OUT], ez[OUT], common_params)
    p7 = point(ex[OUT], ey[OUT], ez[OUT], common_params)
    p8 = point(sx[OUT], ey[OUT], sz[OUT], common_params)


    # PUT PRESSURE ON MIN_W FACE - FACE #1 - SEE FIGURE BELOW
    p1_P = copy.deepcopy(p1)
    p1_P.set_param('PRES', self.Pressure)
    
    p2_P = copy.deepcopy(p2)
    p2_P.set_param('PRES', self.Pressure)

    p3_P = copy.deepcopy(p3)
    p3_P.set_param('PRES', self.Pressure)
    
    p4_P = copy.deepcopy(p4)
    p4_P.set_param('PRES', self.Pressure)

    
    ##                           ^ W               / V
    ##                           |                /       
    ##                8          |             7 /   
    ##               .:........................  
    ##             .' |      F2             .'
    ##         5 .:....................._.-' |
    ##           |    |                 |6   |
    ##           |    |      F5(back)   |    |
    ##           |    |                 |    |  ----> U                 
    ##           | F6 |    F3(front)    | F4 |
    ##           |    |                 |    |
    ##           |    |4                |    |
    ##           |   .'.................|....' 3
    ##           |  /        F1         | .'
    ##           `.:....................|'
    ##           1                      2
    ##                                  
    # FROM CALCULIX MANUAL... (RIGHT HAND NORMALS POINT IN)
    #  face 1: 1-2-3-4 -> FTAGS.MIN_W - 1
    #  face 2: 5-8-7-6 -> FTAGS.MAX_W - 2
    #  face 3: 1-5-6-2 -> FTAGS.MIN_V - 4
    #  face 4: 2-6-7-3 -> FTAGS.MAX_U - 8
    #  face 5: 3-7-8-4 -> FTAGS.MAX_V - 16
    #  face 6: 4-8-5-1 -> FTAGS.MIN_U - 32

    #f1 = face('F1')
    #f1.define_face_from_points(p1, p2, p3, p4)
        
    f1 = face('F1')
    f1.define_face_from_points(p1_P, p2_P, p3_P, p4_P)
        
    f2 = face('F2')
    f2.define_face_from_points(p5, p8, p7, p6)

    f3 = face('F3')
    f3.define_face_from_points(p1, p5, p6, p2)
    
    f4 = face('F4')
    f4.define_face_from_points(p2, p6, p7, p3)
    
    f5 = face('F5')
    f5.define_face_from_points(p3, p7, p8, p4)

    f6 = face('F6')
    f6.define_face_from_points(p4, p8, p5, p1)
    
    bname = 'A'
    blk_a = self.add_block(bname)
    blk_a.set_material(mat)
    blk_a.define_block_from_faces(f1, f2, f3, f4, f5, f6)



#  ---- ---- ---- ---- B BLOCK


    sx = [ hinge_x, hinge_x ]
    ex = [ ring_x, ring_x ]

    sy = [ tys, tys ]
    ey = [ tye, tye ]

    # ALWAYS ADD THE Z THICKNESS TO THE IN SURFACE TO GET THE OUT SURFACE
    sz = [ hinge_z, hinge_z + thref ]
    ez = [ ring_z, ring_z + thref]

    
    p1 = point(sx[IN], sy[IN], sz[IN], common_params)
    p2 = point(ex[IN], sy[IN], ez[IN], common_params)
    p3 = point(ex[IN], ey[IN], ez[IN], common_params)
    p4 = point(sx[IN], ey[IN], sz[IN], common_params)
    p5 = point(sx[OUT], sy[OUT], sz[OUT], common_params)
    p6 = point(ex[OUT], sy[OUT], ez[OUT], common_params)
    p7 = point(ex[OUT], ey[OUT], ez[OUT], common_params)
    p8 = point(sx[OUT], ey[OUT], sz[OUT], common_params)


    # PUT PRESSURE ON MIN_W FACE - FACE #1 - SEE FIGURE BELOW
    p1_P = copy.deepcopy(p1)
    p1_P.set_param('PRES', self.Pressure)
    
    p2_P = copy.deepcopy(p2)
    p2_P.set_param('PRES', self.Pressure)

    p3_P = copy.deepcopy(p3)
    p3_P.set_param('PRES', self.Pressure)
    
    p4_P = copy.deepcopy(p4)
    p4_P.set_param('PRES', self.Pressure)

    
    ##                           ^ W               / V
    ##                           |                /       
    ##                8          |             7 /   
    ##               .:........................  
    ##             .' |      F2             .'
    ##         5 .:....................._.-' |
    ##           |    |                 |6   |
    ##           |    |      F5(back)   |    |
    ##           |    |                 |    |  ----> U                 
    ##           | F6 |    F3(front)    | F4 |
    ##           |    |                 |    |
    ##           |    |4                |    |
    ##           |   .'.................|....' 3
    ##           |  /        F1         | .'
    ##           `.:....................|'
    ##           1                      2
    ##                                  
    # FROM CALCULIX MANUAL... (RIGHT HAND NORMALS POINT IN)
    #  face 1: 1-2-3-4 -> FTAGS.MIN_W - 1
    #  face 2: 5-8-7-6 -> FTAGS.MAX_W - 2
    #  face 3: 1-5-6-2 -> FTAGS.MIN_V - 4
    #  face 4: 2-6-7-3 -> FTAGS.MAX_U - 8
    #  face 5: 3-7-8-4 -> FTAGS.MAX_V - 16
    #  face 6: 4-8-5-1 -> FTAGS.MIN_U - 32

    #f1 = face('F1')
    #f1.define_face_from_points(p1, p2, p3, p4)

    f1 = face('F1')
    f1.define_face_from_points(p1_P, p2_P, p3_P, p4_P)
        
    f2 = face('F2')
    f2.define_face_from_points(p5, p8, p7, p6)

    f3 = face('F3')
    f3.define_face_from_points(p1, p5, p6, p2)
    
    f4 = face('F4')
    f4.define_face_from_points(p2, p6, p7, p3)
    
    f5 = face('F5')
    f5.define_face_from_points(p3, p7, p8, p4)

    f6 = face('F6')
    f6.define_face_from_points(p4, p8, p5, p1)
    
    bname = 'B'
    blk_b = self.add_block(bname)
    blk_b.set_material(mat)
    blk_b.define_block_from_faces(f1, f2, f3, f4, f5, f6)



#  ---- ---- ---- ---- C BLOCK
    
    sx = [ ring_x, ring_x ]
    ex = [ txe, txe ]

    sy = [ tys, tys ]
    ey = [ tye, tye ]

    # ALWAYS ADD THE Z THICKNESS TO THE IN SURFACE TO GET THE OUT SURFACE
    sz = [ ring_z, ring_z + thref ]
    #ez = [ tze + thref/4.0, tze + thref - thref/4.0 ]
    ez = [ tze, tze + thref - thref/2.0 ]

    
    p1 = point(sx[IN], sy[IN], sz[IN], common_params)
    p2 = point(ex[IN], sy[IN], ez[IN], dict( common_params.items() + {'VOL_FRAC':0.6}.items() ) ) # OVERRIDE VOL_FRAC PARAMETER
    p3 = point(ex[IN], ey[IN], ez[IN], dict( common_params.items() + {'VOL_FRAC':0.6}.items() ) )
    p4 = point(sx[IN], ey[IN], sz[IN], common_params)
    p5 = point(sx[OUT], sy[OUT], sz[OUT], common_params)
    p6 = point(ex[OUT], sy[OUT], ez[OUT], dict( common_params.items() + {'VOL_FRAC':0.6}.items() ) )
    p7 = point(ex[OUT], ey[OUT], ez[OUT], dict( common_params.items() + {'VOL_FRAC':0.6}.items() ) )
    p8 = point(sx[OUT], ey[OUT], sz[OUT], common_params)

    # PUT PRESSURE ON MIN_W FACE - FACE #1 - SEE FIGURE BELOW
    p1_P = copy.deepcopy(p1)
    p1_P.set_param('PRES', self.Pressure)
    
    p2_P = copy.deepcopy(p2)
    p2_P.set_param('PRES', self.Pressure)

    p3_P = copy.deepcopy(p3)
    p3_P.set_param('PRES', self.Pressure)
    
    p4_P = copy.deepcopy(p4)
    p4_P.set_param('PRES', self.Pressure)

    
    ##                           ^ W               / V
    ##                           |                /       
    ##                8          |             7 /   
    ##               .:........................  
    ##             .' |      F2             .'
    ##         5 .:....................._.-' |
    ##           |    |                 |6   |
    ##           |    |      F5(back)   |    |
    ##           |    |                 |    |  ----> U                 
    ##           | F6 |    F3(front)    | F4 |
    ##           |    |                 |    |
    ##           |    |4                |    |
    ##           |   .'.................|....' 3
    ##           |  /        F1         | .'
    ##           `.:....................|'
    ##           1                      2
    ##                                  
    # FROM CALCULIX MANUAL... (RIGHT HAND NORMALS POINT IN)
    #  face 1: 1-2-3-4 -> FTAGS.MIN_W - 1
    #  face 2: 5-8-7-6 -> FTAGS.MAX_W - 2
    #  face 3: 1-5-6-2 -> FTAGS.MIN_V - 4
    #  face 4: 2-6-7-3 -> FTAGS.MAX_U - 8
    #  face 5: 3-7-8-4 -> FTAGS.MAX_V - 16
    #  face 6: 4-8-5-1 -> FTAGS.MIN_U - 32

    #f1 = face('F1')
    #f1.define_face_from_points(p1, p2, p3, p4)

    f1 = face('F1')
    f1.define_face_from_points(p1_P, p2_P, p3_P, p4_P)
        
    f2 = face('F2')
    f2.define_face_from_points(p5, p8, p7, p6)
    
    f3 = face('F3')
    f3.define_face_from_points(p1, p5, p6, p2)
    
    f4 = face('F4')
    f4.define_face_from_points(p2, p6, p7, p3)
    
    f5 = face('F5')
    f5.define_face_from_points(p3, p7, p8, p4)

    f6 = face('F6')
    f6.define_face_from_points(p4, p8, p5, p1)
    
    bname = 'C'
    blk_c = self.add_block(bname)
    blk_c.set_material(mat)
    blk_c.define_block_from_faces(f1, f2, f3, f4, f5, f6)



    self.plot_blocks()


    # BLOCK :  +/- 1 parametric CUBE CENTERED AT ZERO
    # EDGE INDICIES (su, sv, sw, eu, ev, ew) ALWAYS POINT FROM LO TO HI VALUE
    # FIRST THREE INDICIES ARE THE START POINT THE LAST THREE ARE THE END POINT INDICIES
    #  i.e.   EDGE(-1, 1, 1, 1, 1, 1) POINTS FROM POINT (-1, 1, 1) TO (1, 1, 1)
    #     THE EDGE VECTOR POINTS DOWN THE U AXIS ALONG THE EDGE WHERE V, ARE W ARE THIER MAXIMUM
    #  ONLY ONE OF THE INDICIES WILL CHANGE FROM (-) TO (+) FROM THE START POINT TO THE END POINT
    #  EACH BLOCK STORES 8 CORNER POINTS AND 12 EDGE VECTORS
    

    #                              NAME,  BLOCK,   GRID_TAG
    #self.add_rigid_join('XRA', blk_x, GTAGS.MIN_W,
    #                         blk_a, GTAGS.MIN_U) 

    #                              BLOCK,   GRID_TAG,    HINGE_AXIS_DIR
    self.add_hinge_join('AHB', blk_a, GTAGS.MAX_U, DIR.V,
                        blk_b, GTAGS.MIN_U, DIR.V) 


#    self.add_butt_join('BBC', blk_b, GTAGS.MAX_U,
#                       blk_c, GTAGS.MIN_U) 

    self.add_butt_join_with_grid_omit_tags('BBC', blk_b, GTAGS.MAX_U, [GTAGS.CENTER_W],
                                           blk_c, GTAGS.MIN_U, [GTAGS.CENTER_W] ) 


#    self.add_rigid_join('BBC', blk_b, GTAGS.MAX_U,
#                       blk_c, GTAGS.MIN_U) 
  
#    self.add_rigid_join_with_grid_omit_tags('BBC', blk_b, GTAGS.MAX_U, [GTAGS.CENTER_W],
#                                            blk_c, GTAGS.MIN_U, [GTAGS.CENTER_W]) 





#    self.add_butt_join('XBA', blk_x, GTAGS.MIN_W,
#                              blk_a, GTAGS.MIN_U) 



    # MUST GENERATE THE MESH BEFORE YOU CAN APPLY THE BOUNDARY CONDITIONS
    self.generate_mesh()

    gla = blk_a.get_grid_list_from_tags(GTAGS.MIN_U)
    gla = blk_a.get_grid_list_from_tags(GTAGS.CENTER_W, gla)
    for g in gla :
        #print g.id
        g.bc.add(DOF.DX)
        g.bc.add(DOF.DZ)
    pass


    
    glb = blk_b.get_grid_list_from_tags(GTAGS.MAX_U)
    glb = blk_b.get_grid_list_from_tags(GTAGS.CENTER_W, glb)
    for g in glb :
        g.bc.add(DOF.DX)
        g.bc.add(DOF.DZ)
    pass

    ## glc = blk_c.get_grid_list_from_tags(GTAGS.MIN_U)
    ## glc = blk_c.get_grid_list_from_tags(GTAGS.CENTER_W, glc)
    ## for g in glc :
    ##     g.bc.add(DOF.DX)
    ##     g.bc.add(DOF.DZ)
    ## pass

    ## glc = blk_c.get_grid_list_from_tags(GTAGS.MAX_U)
    ## for g in glc :
    ##     g.bc.add(DOF.DX)
    ##     g.bc.add(DOF.DZ)
    ## pass


    # AVOID THE FACES THAT HAVE JOINS OR BOUNDARY CONDITIONS
    gcl_a = blk_a.get_grid_list_from_tags(GTAGS.CENTER_V)
    gcl_a = blk_a.get_grid_list_from_not_tags(GTAGS.MAX_U, gcl_a)
    gcl_a = blk_a.get_grid_list_from_not_tags(GTAGS.MIN_U, gcl_a)
    
    gcl_b = blk_b.get_grid_list_from_tags(GTAGS.CENTER_V)
    gcl_b = blk_b.get_grid_list_from_not_tags(GTAGS.MAX_U, gcl_b)
    gcl_b = blk_b.get_grid_list_from_not_tags(GTAGS.MIN_U, gcl_b)
    
    gcl_c = blk_c.get_grid_list_from_tags(GTAGS.CENTER_V)
    gcl_c = blk_c.get_grid_list_from_not_tags(GTAGS.MIN_U, gcl_c)
   
    gcl = gcl_a + gcl_b + gcl_c
    for g in gcl :
        g.bc.add(DOF.DY)
        #g.bc.add(DOF.RX)
    pass


    
    print ' *** M O D E L   B U I L D   D O N E *** '

    #sys.exit(-1)
pass

#####################################################################


anal.analyze(define_model)



