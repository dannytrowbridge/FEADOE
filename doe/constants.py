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
    'DEFINED',
    'EXPLODED'
    ]


BTAGS = namedtuple('BTAGS',
                   block_tags, verbose = False)._make( int(2**i) for i in range(len(block_tags)) )


grid_tags = [
    'CENTER_U',
    'CENTER_V',
    'CENTER_W',
    'MIN_U',
    'MAX_U',
    'MIN_V',
    'MAX_V',
    'MIN_W',
    'MAX_W',
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
    'MERGED',
    'MERGED_FACE',
    'MERGED_EXTRA',
    'REMOVED',
    'EXPANDED_PARENT',
    'EXPANDED',
    'EXPANDED_CHILD',
    'EXPANDED_CLONE',
    'EXPANDED_MIDDLE'
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

## face_tags = [
##     'MAX_U',
##     'MAX_V',
##     'MAX_W',
##     'MIN_U',
##     'MIN_V',
##     'MIN_W'
##     ]
face_tags = [
    'MIN_W',
    'MAX_W',
    'MIN_V',
    'MAX_U',
    'MAX_V',
    'MIN_U',
    'JOIN_FACE_MASTER',
    'JOIN_FACE_SLAVE'
    ]

FTAGS = namedtuple('FTAGS',
                   face_tags, verbose = False)._make( int(2**i) for i in range(len(face_tags)) )



# RIGHT HAND RULE NORMALS POINT IN
# FROM CALCULIX MANUAL...
#  face 1: 1-2-3-4 -> FTAGS.MIN_W - 1
#  face 2: 5-8-7-6 -> FTAGS.MAX_W - 2
#  face 3: 1-5-6-2 -> FTAGS.MIN_V - 4
#  face 4: 2-6-7-3 -> FTAGS.MAX_U - 8
#  face 5: 3-7-8-4 -> FTAGS.MAX_V - 16
#  face 6: 4-8-5-1 -> FTAGS.MIN_U - 32

#face_num = (math.log(etag) / math.log(2.0)) + 1

FACE_GRID_INDICES = [[]] * 6
# WE WILL KEEP THESE VALUES ONES-BASED FOR CLARITY WHEN COMPARING TO THE CALCULIX FACE DEFINITIONS 
FACE_GRID_INDICES[int(math.log(float(FTAGS.MIN_W)) / math.log(2.0))] = [ 1, 2, 3, 4 ] 
FACE_GRID_INDICES[int(math.log(float(FTAGS.MAX_W)) / math.log(2.0))] = [ 5, 8, 7, 6 ] 
FACE_GRID_INDICES[int(math.log(float(FTAGS.MIN_V)) / math.log(2.0))] = [ 1, 5, 6, 2 ] 
FACE_GRID_INDICES[int(math.log(float(FTAGS.MAX_U)) / math.log(2.0))] = [ 2, 6, 7, 3 ] 
FACE_GRID_INDICES[int(math.log(float(FTAGS.MAX_V)) / math.log(2.0))] = [ 3, 7, 8, 4 ] 
FACE_GRID_INDICES[int(math.log(float(FTAGS.MIN_U)) / math.log(2.0))] = [ 4, 8, 5, 1 ] 






join_tags = [
    'MERGE',
    'HINGE'
    ]

JTAGS = namedtuple('JTAGS',
                   join_tags, verbose = False)._make( i+1 for i in range(len(join_tags)) )


