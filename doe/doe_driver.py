import sys
import os
import re
import math
import string
import itertools
from collections import OrderedDict
from ggen_util import print_timing, enum, here, now, interpolate


##########################################################################################################


class indep_var(object) :
    
    def __init__(self, vname, vv = []) :
        self.var = vname
        self.vals = vv
        self.current_index = 0
    pass

        
    #----------------------------------------------------------------------------------------------------

    def make_var_def_string(self) :
        s = self.var + ' = ' + self.vals[current_index]
        return(s)
    pass

pass

##########################################################################################################

class design_of_experiments(object) :

    def __init__ (self, doe_name = '') :
        self.name = doe_name
        self.template = None
        self.indep = OrderedDict()
        self.current_vals = None
        self.combo = None
    pass

    #----------------------------------------------------------------------------------------------------

    def build_combo(self) :
        lol = []
        for k, v in self.indep.items() :
            lol.append(v.vals)
        pass
        self.combo = itertools.product(*lol)
    pass

    #----------------------------------------------------------------------------------------------------

    def add_indep_var(self, ivname, values) :
        indepv = indep_var(ivname, values)
        self.indep[ivname] = indepv
    pass

    #----------------------------------------------------------------------------------------------------

    def get_sample(self) :
        if( self.combo is None ) :
            self.build_combo()
        pass

        try:
            vals = self.combo.next()
            self.current_vals = OrderedDict( zip( self.indep.keys(), vals ) )
        except StopIteration: 
            self.current_vals = None
        pass

        return( self.current_vals )
    pass

    #----------------------------------------------------------------------------------------------------

    def get_var_defs(self) :
        defs = []
        vals = self.get_sample()
        if( vals is None ) : return(defs)
        for k, v in vals.items() :
            s = k + ' = ' + str(v)
            defs.append(s)
        pass
        return(defs)
    pass
        
    
pass

##########################################################################################################
        
if __name__ == '__main__' :

    doe = design_of_experiments()

    

    doe.add_indep_var('r0', [12.0, 13.0, 17.0])
    doe.add_indep_var('r1', [12.1, 13.1, 17.1])
    doe.add_indep_var('r2', [12.2, 13.2, 17.2])
    doe.add_indep_var('r3', [12.3, 13.3, 17.3])

    vals = [0]
    cnt = 0
    while ( len(vals) != 0 ) :
         cnt += 1
         vals = doe.get_var_defs()
         if( len(vals) == 0 ) : break
         print cnt, '  ', vals
         for s in vals :
             exec(s)
         pass
         print 'r0=', r0
         print 'r1=', r1
         print 'r2=', r2
         print 'r3=', r3
    pass


pass
