import sys
import doe.analysis as doe_analysis
import doe.gfea_model as fea
from doe.ggen_util import print_timing, enum, here, now, interpolate
from doe.ggen_sup import point, vector, face
from doe.constants import LOC, DIR, DOF, GTAGS, ETAGS, FTAGS, JTAGS

import math
import copy

anal = doe_analysis.analysis('NOZ_T1')
anal.add_indep_var('Pressure', [-10000.0])
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

    common_params =  { 'PRES': 0.0, 'THICK': 1.0 }

    # FROM CALCULIX MANUAL...
    #  face 1: 1-2-3-4 -> FTAGS.MIN_W - 1
    #  face 2: 5-8-7-6 -> FTAGS.MAX_W - 2
    #  face 3: 1-5-6-2 -> FTAGS.MIN_V - 4
    #  face 4: 2-6-7-3 -> FTAGS.MAX_U - 8
    #  face 5: 3-7-8-4 -> FTAGS.MAX_V - 16
    #  face 6: 4-8-5-1 -> FTAGS.MIN_U - 32

    # Must use faces (not points) so we can change pressures, etc. independently on each face  

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
## FACE NORMALS POINT IN


    ## xlox = 0.0
    ## xhix = 1.0
    ## xoffx = -0.5
    
    ## xloy = -0.5
    ## xhiy = 0.5
    ## xoffy = 0.0
    
    ## xloz = 0.0
    ## xhiz = 4.0
    ## xoffz = 0.0
    
    ## px1 = point(xlox + xoffx, xloy + xoffy, xloz + xoffz, common_params)
    ## px2 = point(xhix + xoffx, xloy + xoffy, xloz + xoffz, common_params)
    ## px3 = point(xhix + xoffx, xhiy + xoffy, xloz + xoffz, common_params)
    ## px4 = point(xlox + xoffx, xhiy + xoffy, xloz + xoffz, common_params)
    ## px5 = point(xlox + xoffx, xloy + xoffy, xhiz + xoffz, common_params)
    ## px6 = point(xhix + xoffx, xloy + xoffy, xhiz + xoffz, common_params)
    ## px7 = point(xhix + xoffx, xhiy + xoffy, xhiz + xoffz, common_params)
    ## px8 = point(xlox + xoffx, xhiy + xoffy, xhiz + xoffz, common_params)


    ## bname = 'X'
    ## blk_x = self.add_block(bname)
    ## blk_x.set_material(mat)
    ## #blk_x.force_element_count([2, 12, 12])
    ## # HAVE TO CALL DEFINE BLOCK ROUTINE LAST 
    ## blk_x.define_block_from_points(px1, px2, px3, px4, px5, px6, px7, px8)

    txs = 0.0
    txe = 21.0

    tys = -2.0
    tye = 2.0

    tzs = 10.0
    tze = 11.0
    thref = 1.0

    hinge_x = 9.0
    hinge_z = 5.0

    ring_x = 15.0
    ring_z = interpolate(hinge_z, hinge_x, tze, txe, ring_x)

       
#  ---- ---- ---- ---- A BLOCK

    sx = [ txs, txs ]
    ex = [ hinge_x, hinge_x ]

    sy = [ tys, tys ]
    ey = [ tye, tye ]

    # ALWAYS ADD THE Z THICKNESS TO THE OUT SURFACE
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


    # PUT PRESSURE ON MAX_W FACE - FACE #2 - SEE FIGURE BELOW
    p5_P = copy.deepcopy(p5)
    p5_P.set_param('PRES', self.Pressure)
    
    p8_P = copy.deepcopy(p8)
    p8_P.set_param('PRES', self.Pressure)

    p7_P = copy.deepcopy(p7)
    p7_P.set_param('PRES', self.Pressure)
    
    p6_P = copy.deepcopy(p6)
    p6_P.set_param('PRES', self.Pressure)

    
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

    f1 = face('F1')
    f1.define_face_from_points(p1, p2, p3, p4)
        
    #f2 = face('F2')
    #f2.define_face_from_points(p5, p8, p7, p6)

    f2 = face('F2')
    f2.define_face_from_points(p5_P, p8_P, p7_P, p6_P)

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


    lox = hinge_x
    hix = ring_x

    loy = tys
    hiy = tye

    
    loz = hinge_z
    hiz = ring_z

    
    p1 = point( lox, loy, loz + 0.5, { 'PRES': self.Pressure, 'THICK': 1.0 })
    p2 = point( hix, loy, hiz + 0.4, { 'PRES': self.Pressure * 2.0, 'THICK': 0.8 })
    p3 = point( hix, hiy, hiz + 0.4, { 'PRES': self.Pressure * 2.0, 'THICK': 0.8 })
    p4 = point( lox, hiy, loz + 0.5, { 'PRES': self.Pressure, 'THICK': 1.0 })


    bname = 'B'
    blk_b = self.add_block(bname)
    blk_b.set_material(mat)
    blk_b.define_block_from_mid_plane_points(p1, p2, p3, p4)


#  ---- ---- ---- ---- C BLOCK
    
    sx = [ ring_x, ring_x ]
    ex = [ txe, txe ]

    sy = [ tys, tys ]
    ey = [ tye, tye ]

    # ALWAYS ADD THE Z THICKNESS TO THE OUT SURFACE
    sz = [ ring_z, ring_z + 0.8 ]
    ez = [ tze, tze + 0.8 ]

    
    p1 = point(sx[IN], sy[IN], sz[IN], common_params)
    p2 = point(ex[IN], sy[IN], ez[IN], common_params)
    p3 = point(ex[IN], ey[IN], ez[IN], common_params)
    p4 = point(sx[IN], ey[IN], sz[IN], common_params)
    p5 = point(sx[OUT], sy[OUT], sz[OUT], common_params)
    p6 = point(ex[OUT], sy[OUT], ez[OUT], common_params)
    p7 = point(ex[OUT], ey[OUT], ez[OUT], common_params)
    p8 = point(sx[OUT], ey[OUT], sz[OUT], common_params)

    # PUT PRESSURE ON MAX_W FACE - FACE #2 - SEE FIGURE BELOW
    p5_P = copy.deepcopy(p5)
    p5_P.set_param('PRES', self.Pressure * 2.0)
    
    p8_P = copy.deepcopy(p8)
    p8_P.set_param('PRES', self.Pressure * 2.0)

    p7_P = copy.deepcopy(p7)
    p7_P.set_param('PRES', self.Pressure)
    
    p6_P = copy.deepcopy(p6)
    p6_P.set_param('PRES', self.Pressure)

    
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

    f1 = face('F1')
    f1.define_face_from_points(p1, p2, p3, p4)
        
    #f2 = face('F2')
    #f2.define_face_from_points(p5, p8, p7, p6)

    f2 = face('F2')
    f2.define_face_from_points(p5_P, p8_P, p7_P, p6_P)

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


# WE CAN'T USE BUTT HERE BECAUSE THE MAX_U FACE OF BLOCK B IS NOT COPLANER WITH THE
# MIN_U FACE OF BLOCK C
#    self.add_butt_join('BBC', blk_b, GTAGS.MAX_U,
#                       blk_c, GTAGS.MIN_U) 

    self.add_rigid_join('BBC', blk_b, GTAGS.MAX_U,
                       blk_c, GTAGS.MIN_U) 





#    self.add_butt_join('XBA', blk_x, GTAGS.MIN_W,
#                              blk_a, GTAGS.MIN_U) 



    # MUST GENERATE THE MESH BEFORE YOU CAN APPLY THE BOUNDARY CONDITIONS
    self.generate_mesh()

    gla = blk_a.get_grid_list_from_tags(GTAGS.MIN_U)
    gla = blk_a.get_grid_list_from_tags(GTAGS.CENTER_W, gla)
    #print 'FACE POINTS...'
    for g in gla :
        #print g.id
        g.bc.add(DOF.DX)
        g.bc.add(DOF.DZ)
    pass


    
    glb = blk_b.get_grid_list_from_tags(GTAGS.MAX_U)
    for g in glb :
        g.bc.add(DOF.DX)
        g.bc.add(DOF.DZ)
    pass

    ## glc = blk_c.get_grid_list_from_tags(GTAGS.MIN_U)
    ## for g in glc :
    ##     g.bc.add(DOF.DX)
    ##     g.bc.add(DOF.DZ)
    ## pass


    ## gcl_a = blk_a.get_grid_list_from_tags(GTAGS.CENTER_V)
    ## gcl_b = blk_b.get_grid_list_from_tags(GTAGS.CENTER_V)
    ## gcl_c = blk_c.get_grid_list_from_tags(GTAGS.CENTER_V)
    ## gcl = gcl_a + gcl_b + gcl_c
    ## for g in gcl :
    ##     g.bc.add(DOF.DY)
    ##     #g.bc.add(DOF.RX)
    ## pass


    
    print ' *** M O D E L   B U I L D   D O N E *** '

    #sys.exit(-1)
pass

#####################################################################


anal.analyze(define_model)



