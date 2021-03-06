import math
from collections import namedtuple
from ggen_util import enum

MATS_FROM_NPSS = 'MATS_FROM_NPSS_RAW.xml'

BIG_REAL = 1.0E120
DEFAULT_INITIAL_TEMP = 70.0
DEFAULT_REFERENCE_TEMP = 0.0
TOL = 1.0E-8

MIN = 0
AVG = 1
MAX = 2

mat_prop_tags = [
    'DENSITY',
    'YM11',
    'YM22',
    'YM33',
    'PR12',
    'PR13',
    'PR23',
    'SM12',
    'SM13',
    'SM23',
    'TEC11',
    'TEC22',
    'TEC33',
    'SY',
    'SU'
    ]

#MAT = enum(mat_prop_tags);
#
#PROPS = namedtuple('PROPS',
#                   mat_prop_tags, verbose = False)._make( i for i in range(len(mat_prop_tags)) )


loc_tags = [
    'X',
    'Y',
    'Z'
    ]

LOC = enum(loc_tags);

dir_tags = [
    'U',
    'V',
    'W'
    ]

DIR = enum(dir_tags);


dof_tags = [
    'DX',
    'DY',
    'DZ',
    'RX',
    'RY',
    'RZ'
    ]

DOF = namedtuple('DOF',
                   dof_tags, verbose = False)._make( i+1 for i in range(len(dof_tags)) )


block_tags = [
    'POINTS',
    'FACES',
    'EXPLODED'
    ]


BTAGS = namedtuple('BTAGS',
                   block_tags, verbose = False)._make( int(2**i) for i in range(len(block_tags)) )

# THE FIRST SIX TERMS OF grid_tags AND face_tags ARE THE SAME FOR CONVIENCE

grid_tags = [
    'MIN_W',
    'MAX_W',
    'MIN_V',
    'MAX_U',
    'MAX_V',
    'MIN_U',
    'CENTER_U',
    'CENTER_V',
    'CENTER_W',
    'ANCILLARY',
    'PENDING_DELETE',
    'ELEMENT_MID_SIDE',
    'ELEMENT_MID_SIDE_U',
    'ELEMENT_MID_SIDE_V',
    'ELEMENT_MID_SIDE_W',
    'ELEMENT_MID_EDGE',
    'ELEMENT_MID_EDGE_U',
    'ELEMENT_MID_EDGE_V',
    'ELEMENT_MID_EDGE_W',
    'ELEMENT_MID_FACE',
    'ELEMENT_MID_FACE_UV',
    'ELEMENT_MID_FACE_VW',
    'ELEMENT_MID_FACE_UW',
    'ELEMENT_CORNER',
    'ELEMENT_MID_VOL',
    'JOIN_MASTER',
    'JOIN_SLAVE',
    'HINGED',
    'HINGED_FACE',
    'HINGED_EXTRA',
    'BUTT',
    'BUTT_FACE',
    'BUTT_EXTRA',
    'RIGID',
    'RIGID_FACE',
    'RIGID_EXTRA',
    'REMOVED'
    ]

GTAGS = namedtuple('GTAGS',
                   grid_tags, verbose = False)._make( int(2**i) for i in range(len(grid_tags)) )


edge_tags = [
    'MAX_U',
    'MAX_V',
    'MAX_W',
    'MIN_U',
    'MIN_V',
    'MIN_W',
    'MAX_U_MAX_V',
    'MAX_U_MIN_V',
    'MIN_U_MAX_V',
    'MIN_U_MIN_V',
    'MAX_V_MAX_W',
    'MAX_V_MIN_W',
    'MIN_V_MAX_W',
    'MIN_V_MIN_W',
    'MAX_U_MAX_W',
    'MAX_U_MIN_W',
    'MIN_U_MAX_W',
    'MIN_U_MIN_W'
    ]

ETAGS = namedtuple('ETAGS',
                   edge_tags, verbose = False)._make( int(2**i) for i in range(len(edge_tags)) )

face_tags = [
    'MIN_W',
    'MAX_W',
    'MIN_V',
    'MAX_U',
    'MAX_V',
    'MIN_U',
    'MID_U',
    'MID_V',
    'MID_W',
    'JOIN_FACE_MASTER',
    'JOIN_FACE_SLAVE'
    ]

FTAGS = namedtuple('FTAGS',
                   face_tags, verbose = False)._make( int(2**i) for i in range(len(face_tags)) )


# ONES BASED FACE INDEX NUMBER FROM CALCULIX
#FACE = [ 'DUMMY', FTAGS.MIN_W, FTAGS.MAX_W, FTAGS.MIN_V, FTAGS.MAX_U, FTAGS.MAX_V, FTAGS.MIN_U ]
def FACE_NUM(ftag) : return(int((math.log(ftag) / math.log(2.0)) + 1))
## FACE_GRID_INDICES[int(math.log(float(FTAGS.MIN_W)) / math.log(2.0))] = [ 1, 2, 3, 4 ]
# ONE'S BASED SO IT MATCHES CALCULIX 
CORNER_POINT_INDEX = [
    ( 0,  0,  0),   # DUMMY
    (-1, -1, -1),   # 1
    ( 1, -1, -1),   # 2
    ( 1,  1, -1),   # 3
    (-1,  1, -1),   # 4
    (-1, -1,  1),   # 5
    ( 1, -1,  1),   # 6
    ( 1,  1,  1),   # 7
    (-1,  1,  1)    # 8
    ]


# RIGHT HAND RULE NORMALS POINT IN
# FROM CALCULIX MANUAL...
#  face 1: 1-2-3-4 -> FTAGS.MIN_W - 1
#  face 2: 5-8-7-6 -> FTAGS.MAX_W - 2
#  face 3: 1-5-6-2 -> FTAGS.MIN_V - 4
#  face 4: 2-6-7-3 -> FTAGS.MAX_U - 8
#  face 5: 3-7-8-4 -> FTAGS.MAX_V - 16
#  face 6: 4-8-5-1 -> FTAGS.MIN_U - 32


FACE_GRID_INDICES = [[]] * 6
# WE WILL KEEP THESE >VALUES< ONES-BASED FOR CLARITY WHEN COMPARING TO THE CALCULIX FACE DEFINITIONS 
## FACE_GRID_INDICES[int(math.log(float(FTAGS.MIN_W)) / math.log(2.0))] = [ 1, 2, 3, 4 ] 
## FACE_GRID_INDICES[int(math.log(float(FTAGS.MAX_W)) / math.log(2.0))] = [ 5, 8, 7, 6 ] 
## FACE_GRID_INDICES[int(math.log(float(FTAGS.MIN_V)) / math.log(2.0))] = [ 1, 5, 6, 2 ] 
## FACE_GRID_INDICES[int(math.log(float(FTAGS.MAX_U)) / math.log(2.0))] = [ 2, 6, 7, 3 ] 
## FACE_GRID_INDICES[int(math.log(float(FTAGS.MAX_V)) / math.log(2.0))] = [ 3, 7, 8, 4 ] 
## FACE_GRID_INDICES[int(math.log(float(FTAGS.MIN_U)) / math.log(2.0))] = [ 4, 8, 5, 1 ] 

FACE_GRID_INDICES[ FACE_NUM(FTAGS.MIN_W) - 1 ] = [ 1, 2, 3, 4 ] 
FACE_GRID_INDICES[ FACE_NUM(FTAGS.MAX_W) - 1 ] = [ 5, 8, 7, 6 ] 
FACE_GRID_INDICES[ FACE_NUM(FTAGS.MIN_V) - 1 ] = [ 1, 5, 6, 2 ] 
FACE_GRID_INDICES[ FACE_NUM(FTAGS.MAX_U) - 1 ] = [ 2, 6, 7, 3 ] 
FACE_GRID_INDICES[ FACE_NUM(FTAGS.MAX_V) - 1 ] = [ 3, 7, 8, 4 ] 
FACE_GRID_INDICES[ FACE_NUM(FTAGS.MIN_U) - 1 ] = [ 4, 8, 5, 1 ] 


join_tags = [
    'BUTT',
    'RIGID',
    'HINGE'
    ]

JTAGS = namedtuple('JTAGS',
                   join_tags, verbose = False)._make( i+1 for i in range(len(join_tags)) )


