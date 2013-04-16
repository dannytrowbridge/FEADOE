import sys
import os
import re
import math
import glob
import zipfile
import types
import string
import itertools
from collections import OrderedDict
from WATE_waiter import echo_and_execute_cmd
import gfea_model as fea
from doe_driver import design_of_experiments
from ggen_util import print_timing, enum, here, now, interpolate

SOLVER = 'E:/DEV/pywork/FEADOE/calculix/do_ccx.bat'

if( 'PDRIVE' in os.environ ) :
    PDRIVE = os.environ['PDRIVE']
    #SOLVER = PDRIVE + '/DEV/pywork/FEADOE/calculix/do_ccx.bat'
    SOLVER = PDRIVE + '/DEV/pywork/FEADOE/calculix/do_ccx.bat'
pass

if( 'FEADOE_SOLVER' in os.environ ) :
    SOLVER = os.environ['FEADOE_SOLVER']
pass

#print 'THIS MODULE:', __name__, ' LIVES AT', os.path.realpath(__file__)

##########################################################################################################

class analysis(object) :
    
    #----------------------------------------------------------------------------------------------------

    def __init__(self, name) :
        self.name = name
        self.model = None
        self_model_def_function = None
        self.doe = design_of_experiments()
        self.solver_root_name = ''
        self.solver_in_name = ''
        self.solver_out_name = ''
        self.solver_out_frd_name = ''
        self.solver_grid_map_name = ''
        self.pickle_file = ''
        self.analysis_count = 0
        self.min_results = {}
        self.max_results = {}
        self.solver_name = SOLVER

        #self.doe_sum = OrderedDict()

        self.zip_name = self.name + '_MODELS.zip'
        self.summary_name = self.name + '.csv'
        #self.doe_summary_name = self.name + '_DOE_SUMMARY.txt'
        self.user_name = 'UNDEFINED'
        if( os.environ.has_key("USERNAME") ) :
            self.user_name = os.environ['USERNAME']
        pass        

        self.results = OrderedDict()

    pass

    #----------------------------------------------------------------------------------------------------
    # http://stackoverflow.com/questions/2484215/how-do-i-override-a-parent-classs-functions-in-python
    # http://stackoverflow.com/questions/972/adding-a-method-to-an-existing-object
    def set_model_def(self) :
        print '>>> TOP >>> ', here()
        #self.model.build_model = new.instancemethod(md_func, self.model, self.model.__class__)
        self.model.build_model = types.MethodType(self.model_def_function, self.model)
    pass

    #----------------------------------------------------------------------------------------------------

    def add_indep_var(self, ivname, values) :
        print '>>> TOP >>> ', here()
        self.doe.add_indep_var(ivname, values)
    pass

    #----------------------------------------------------------------------------------------------------
    # UPDATES THE DOE VARIABLES TO THE NEXT SET OF VALUES - 'None' MEANS THE DOE IS DONE
    def update_vars(self) :
        rc = 0
        current_vals = self.doe.get_sample()
        if( current_vals is None ) :
            rc = -1
        pass
        print 'CURRENT VALS = ', current_vals
        #self.doe_sum[self.analysis_count] = current_vals       
        return(rc)
    pass

    #----------------------------------------------------------------------------------------------------

    def insert_vars_into_model(self) :
        self.model.header_text_lines = []
        for k, v in self.doe.current_vals.items() :
            setattr(self.model, k, v)
            s = '  ' + k + ' = ' + str(v) 
            self.model.header_text_lines.append(s)
        pass
    pass

    #----------------------------------------------------------------------------------------------------

    def analyze(self, md_func) :
        self.model_def_function = md_func
        self.analysis_count = 1
        rc = self.update_vars()
        while( rc == 0 ) :
            print 'ANALYSIS COUNT = ', self.analysis_count

            #self.model.__init__()
            del self.model
            self.model = fea.model()
            self.set_model_def()
            
            self.insert_vars_into_model()
            self.define_model()
            self.submit_to_solver()
            self.extract_results()
            rc = self.update_vars()
            if( rc == 0 ) :
                self.analysis_count += 1
            pass
        pass
    
        self.print_summary()
        self.zip_things_up()
    pass

    #----------------------------------------------------------------------------------------------------

    def define_model(self) :
        self.model.build_model()
    pass

    #----------------------------------------------------------------------------------------------------

    def submit_to_solver(self) :
        self.solver_root_name = self.name + '_X' + str(self.analysis_count)
        self.solver_in_name = self.solver_root_name + '.inp'
        self.solver_out_name = self.solver_root_name + '.dat'
        self.solver_out_frd_name = self.solver_root_name + '.frd'
        self.solver_grid_map_name = self.solver_root_name + '.XML'
        self.pickle_file = self.solver_root_name + 'pickle'
        #self.model.save_abq_QUAD(self.solver_in_name)
        self.model.save_abq(self.solver_in_name)
        self.solver_command = self.solver_name + ' ' + self.solver_root_name
        echo_and_execute_cmd(self.solver_command) # THIS WAITS
    pass

    #----------------------------------------------------------------------------------------------------

    def extract_results(self) :
        model_3d_name = self.solver_root_name + '_3D.inp'
        #self.model.read_stress_results_from_dot_frd_QUAD(self.solver_out_frd_name,
        #                                            self.solver_grid_map_name,
        #                                            model_3d_name)
        self.model.read_stress_results_from_dot_frd(self.solver_out_frd_name)



        #self.model.read_stress_results_from_dot_dat(self.solver_out_name)
        # NEXT LINE THROWS... TypeError: can't pickle instancemethod objects
        # NG del self.model.build_model
        #self.model.save(self.pickle_file)
        
        extra = {}
        extra['AID'] = self.analysis_count;
        extra['ANALYSIS_NAME'] = self.solver_root_name
        for k, v in self.doe.current_vals.items() :
            extra[k] = v
        pass
    
        # self.results is an OrderedDict
        # self.model.max_results is a dict
        # self.model.avg_results is a dict
        # self.model.min_results is a dict
        # extra is a dict
        #                                              0                      1                        2               3
        self.results[self.analysis_count] = [ self.model.max_results, self.model.avg_results, self.model.min_results, extra ]  

    
        for k, v in self.model.max_results.items() :
            if( k not in self.max_results ) :
                 self.max_results[k] = [v, extra]
            pass
            #print 'k = ', k, 'CURRENT MAX = ', self.max_results[k][0], '  LOOKING AT:', v
            if( self.max_results[k][0] < v ) :
                #print '    NEW MAX =', v
                self.max_results[k] = [v, extra]
            pass
        pass
    
        for k, v in self.model.min_results.items() :
            if( k not in self.min_results ) :
                 self.min_results[k] = [v, extra]
            pass
            #print 'k = ', k, 'CURRENT MIN = ', self.min_results[k][0], '  LOOKING AT:', v
            if( self.min_results[k][0] > v ) :
                #print '    NEW MIN =', v
                self.min_results[k] = [v, extra]
            pass
        pass
                
    pass

    #----------------------------------------------------------------------------------------------------

    def write_csv_individual_analysis_vals(self, fp, index, tag) :

        fp.write(tag + ',')
        #print '\nTAG var = ', tag
        for k, v in self.results.items() :
            #print 'k=', k
            #print 'v=', v
            rrs = v[index]
            s = ','
            #print 'rrs =', rrs
            if tag in rrs :
                #print 'rrs[tag] =', rrs[tag]
                s = str(rrs[tag]) + ','
            else :
                s = 'NA,'
            pass
            #print 's=', s
            fp.write(s)
        pass
        fp.write('\n')
    pass


    #----------------------------------------------------------------------------------------------------

    def write_csv_overal_doe_vals(self, fp, index, stag, results_dict, tag) :

        fp.write(tag + ',')
        #print '\nTAG var = ', tag
        v = results_dict[tag]
        fp.write( stag + '_' + tag + ', ' + str(v[0]))
        #print 'VAR=', tag, '   VAL=', v[0]
        
        # EXTRA DATA CORRESPONDING TO THE ANALYSIS CASE
        for kk, vv in v[1].items() :
            fp.write( ', ' + kk + ', ' + str(vv))
            #print '     ', kk, '  =  ', vv,
        pass

        fp.write('\n')
    pass


    #----------------------------------------------------------------------------------------------------

    def print_summary(self) :

        MAX = 0
        AVG = 1
        MIN = 2
        EXTRA = 3
        
        id_re = re.compile(r'_ID[ ]*$', re.IGNORECASE )

        b2b_brace_re = re.compile( r'[}][ ]*[{]', re.IGNORECASE )

        infile = sys.modules['__main__'].__file__

        fp = open(self.summary_name, 'w')
        fp.write('#, DOE RESULTS FROM ANLAYSIS "' + self.name + '"\n')
        fp.write('#, THE ANALYSIS TOOK PLACE ON: ' + now() + '\n')
        fp.write('#, THE INPUT FILE WAS: ' + infile + '\n')
        fp.write('#, PERFORMED BY: ' + self.user_name + '\n')
        fp.write('#, ANALYSIS FILES ARE STORED IN : ' + self.zip_name + '\n')
        fp.write('#\n#, INDEPENDENT VARIABLES:' + '\n')

        for k, iv in self.doe.indep.items() :
            f = ' {}' * len(iv.vals)
            f = b2b_brace_re.sub(r'},{', f.strip())
            f = f.format(*iv.vals)
            fp.write('#,' + k + ',' + f + '\n')
            #fp.write('#   ' + k + ' = ' + str(iv.vals) + '\n')
        pass

        fp.write('\n')
        fp.write('#'*80 +'\n')
        fp.write('#, THE NUMBER OF ANALYSES IN DOE = ' + str(self.analysis_count) + '\n')


        fp.write('#'*80 +'\n')

        kk = self.max_results.keys()
        kk.sort()
        #print 'kk=',kk
        #print 'MAX VALUES:'
        for k in kk :
            # SKIP IDS - HANDLED EXPLICTLY BELOW SO THEY APPEAR TOGETHER BELOW EACH LINE OF CORRESPONDING DATA
            if( id_re.search(k) >= 0 ) : continue
            
            self.write_csv_overal_doe_vals(fp, MAX, 'MAX', self.max_results, k)

            kid = k + '_ID'
            if kid in kk :
                self.write_csv_overal_doe_vals(fp, MAX, 'MAX', self.max_results, kid)
            pass      


##             v = self.max_results[k]
##             print 'k=', k
##             print 'v=', v
##             fp.write( 'MAX_' + k + ', ' + str(v[0]))
##             print 'VAR=', k, '   VAL=', v[0]    
##             for kk, vv in v[1].items() :
##                 fp.write( ', ' + kk + ', ' + str(vv))
##                 print '     ', kk, '  =  ', vv,
##             pass
##             fp.write('\n')
##             print
        pass

        fp.write('='*80 + '\n')

        kk = self.min_results.keys()
        kk.sort()
        #print 'MIN VALUES:'
        for k in kk :
            
            # SKIP IDS - HANDLED EXPLICTLY BELOW SO THEY APPEAR TOGETHER BELOW EACH LINE OF CORRESPONDING DATA
            if( id_re.search(k) >= 0 ) : continue
            
            self.write_csv_overal_doe_vals(fp, MIN, 'MIN', self.min_results, k)

            kid = k + '_ID'
            if kid in kk :
                self.write_csv_overal_doe_vals(fp, MIN, 'MIN', self.min_results, kid)
            pass      

##             v = self.min_results[k]
##             fp.write( 'MIN_' + k + ', ' + str(v[0]))
##             print 'VAR=', k, '   VAL=', v[0]
##             for kk, vv in v[1].items() :
##                 fp.write( ', ' + kk + ', ' + str(vv))
##                 print '     ', kk, '  =  ', vv,
##             pass
##             fp.write('\n')
##             print

        pass
    
        fp.write('#'*80 +'\n')
        fp.write('\n\n')

        ex_tags = set()
        res_tags = set()


        
        # *********                      0                       1                     2                 3
        # self.results ORDER [ self.model.max_results, self.model.avg_results, self.model.min_results, extra ]  
        
        for k, v in self.results.items() :
            for i in range(3) :
                rr = v[i]
                for rk, rv in rr.items() : 
                    res_tags.add(rk) # res_tags IS A SET() SO NO DUPLICATES WILL BE ALLOWED
                pass
            pass

            ex = v[EXTRA] # extra
            for ek, ev in ex.items() : 
                ex_tags.add(ek)
            pass

        pass

        res_tags = list(res_tags)
        res_tags.sort()
        ex_tags = list(ex_tags)
        ex_tags.sort()

        #print
        #print
        #print '\nINDEPENDENT VARIABLES:\n'
        

        fp.write('\nINDEPENDENT VARIABLES:\n')
        for t in ex_tags :
            fp.write(t + ',')
            #print '\nT var = ', t
            for k, v in self.results.items() :
                #print 'k=', k
                #print 'v=', v
                
                rrs = v[EXTRA]
                s = ','
                #print 'rrs =', rrs
                if t in rrs :
                    #print 't = ', t
                    #print 'rrs[t] =', rrs[t]
                    s = str(rrs[t]) + ','
                pass
                #print 's=', s
                fp.write(s)
            pass
            fp.write('\n')
        pass


        #print
        #print
        #print '\nMAX RESULTS:\n'





        
        fp.write('\n')

        fp.write('MAX RESULTS:\n')
        for t in res_tags :
            # SKIP IDS - HANDLED EXPLICTLY BELOW SO THEY APPEAR TOGETHER BELOW THE DATA
            if( id_re.search(t) >= 0 ) : continue
            
            self.write_csv_individual_analysis_vals(fp, MAX, t)

            tt = t + '_ID'
            if tt in res_tags :
                self.write_csv_individual_analysis_vals(fp, MAX, tt)
            pass      
                
        pass

        fp.write('\n')

        fp.write('MIN RESULTS:\n')
        for t in res_tags :

            # SKIP IDS - HANDLED EXPLICTLY BELOW SO THEY APPEAR TOGETHER BELOW THE DATA
            if( id_re.search(t) >= 0 ) : continue
            
            self.write_csv_individual_analysis_vals(fp, MIN, t)

            tt = t + '_ID'
            if tt in res_tags :
                self.write_csv_individual_analysis_vals(fp, MIN, tt)
            pass      

        pass

        fp.write('\n')

        fp.write('AVERAGE RESULTS:\n')
        for t in res_tags :
            # SKIP IDS - HANDLED EXPLICTLY BELOW SO THEY APPEAR TOGETHER BELOW THE DATA
            if( id_re.search(t) >= 0 ) : continue
            
            self.write_csv_individual_analysis_vals(fp, AVG, t)

            tt = t + '_ID'
            if tt in res_tags :
                self.write_csv_individual_analysis_vals(fp, AVG, tt)
            pass      
        pass

    
        fp.close()
    pass
    
    #----------------------------------------------------------------------------------------------------

    def zip_things_up(self) :
        zz = zipfile.ZipFile(self.zip_name, 'w', zipfile.ZIP_DEFLATED, True)
        fl = glob.glob(self.name + '*')
        for f in fl :
            if( f == self.zip_name ) : continue
            print "ZIPPING: " + f
            zz.write(f)
        pass
        zz.close()

        fl = glob.glob(self.name + '_X*')
        for f in fl :
            if( f == self.zip_name ) : continue
            try:
                #os.remove(f)
                print "TESTING ON HOLD: DELETING: " + f
            except OSError, e:
                if e.errno <> errno.ENOENT :
                    raise
                pass
            pass
        pass
    pass

    #----------------------------------------------------------------------------------------------------

pass

##########################################################################################################
