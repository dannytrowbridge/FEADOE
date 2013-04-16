import sys
import os
import time
from datetime import date

##################################################################################################/
# TIMING DECORATOR

def print_timing(func):
    def wrapper(*arg):
        t1 = time.time()
        res = func(*arg)
        t2 = time.time()
        print 'TIMING >>> %s took %0.3f ms' % (func.func_name, (t2-t1)*1000.0)
        return res
    return wrapper
pass
                            
##################################################################################################/
## http://code.activestate.com/recipes/413486/
## http://stackoverflow.com/questions/36932/whats-the-best-way-to-implement-an-enum-in-python
##
## The single asterisk form (*args) is used to pass a non-keyworded, variable-length argument list,
## and the double asterisk form is used to pass a keyworded, variable-length argument list.
## http://www.saltycrane.com/blog/2008/01/how-to-use-args-and-kwargs-in-python/
##
## http://stackoverflow.com/questions/3394835/args-and-kwargs
##
## type(name, bases, dict)
##   Return a new type object. This is essentially a dynamic form of the class statement.
##   The name string is the class name and becomes the __name__ attribute;
##   the bases tuple itemizes the base classes and becomes the __bases__ attribute;
##   and the dict dictionary is the namespace containing definitions
##   for class body and becomes the __dict__ attribute. For example,
##   the following two statements create identical type objects:
##>>> class X(object):
##...     a = 1
##...
##>>> X = type('X', (object,), dict(a=1))
## http://docs.python.org/library/functions.html#type

def enum(vlist, *sequential, **named):
    v = len(vlist)
    enums = {}
    if( v > 0 ) :
        enums = dict(zip(vlist, range(len(vlist))) )
    pass

    #Update the dictionary with the key/value pairs from other, overwriting existing keys
    enums.update(dict(zip(sequential, range(v, v+len(sequential))), **named))

    ## http://docs.python.org/library/functions.html#type
    ee = type('Enum', (), enums)
    # ADD THE DICTIONARY TO THE NEW CLASS - THIS IS A LITTLE REDUNDANT
    ee.vals = enums
    return ee
pass

##################################################################################################/
## http://code.activestate.com/recipes/52315/

def here(v = []):
    s = ""
    try:
        fr = sys._getframe(1)
        s = fr.f_code.co_name + ' L:' + str(fr.f_lineno) + ' ' + fr.f_code.co_filename + ' '
        #s = s + ' |' + str(fr.f_builtins["__self__"]) + '| '
        #s = s + ' |' + fr.f_builtins + '| '
        #s = s + ' |' + fr.f_builtins.__name__ + '| '
        if( len(v) != 0 ) :
            ss = "VARS:: ";
            for var in v:
                val = str(fr.f_locals[var])
                ss = ss + ' ' + var + '=' + val
            pass
            s = s + ss
        pass    
    except Exception as e :
        pass
    pass
    return s
pass

##################################################################################################/
# TBD - Replace old style string formatting

def now() :
  endings = ['st', 'nd', 'rd'] + 17 * ['th'] \
        + ['st', 'nd', 'rd'] + 7 * ['th'] \
        + ['st']

  now = date.today()

  day_of_month = int(now.strftime("%d"))
  now_day_of_week_name = now.strftime("%A")
  now_month_of_the_year_name = now.strftime("%B")
  now_year = now.strftime("%Y")

  now_date = "%s, %s %d%s, %s" % ( now_day_of_week_name, now_month_of_the_year_name, day_of_month, endings[day_of_month-1], now_year)

  here_time = time.localtime()
  fmt = '%I:%M:%S %p'
  now_here_time = time.strftime(fmt, here_time)
  ss = "%s at %s" % (now_date, now_here_time)
  return(ss)
pass

##################################################################################################/

#  RETURNS v CORRESPONDING TO t
def interpolate( lowv, lowt, uppv, uppt, t) :
    return( ( ( uppv - lowv ) / ( uppt - lowt ) ) * ( t - lowt ) + lowv );
pass

##################################################################################################/
