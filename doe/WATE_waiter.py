# ====================================================================
import thread
import time
import os
import sys
#DEPRECATED IN 2.6 > import sets
#import timeit

from datetime import date

# NOTE WE ARE NOT USING THE MUTEX PYTHON OBJECT
#  IT PROBABLY WOULD HAVE BEEN BETTER IF WE DID
#  BUT IN ESSENCE WE ARE DOING THE SAME THING IT DOES.

# ====================================================================

def print_timing(func):
    def wrapper(*arg):
        t1 = time.time()
        res = func(*arg)
        t2 = time.time()
        print '%s took %0.3f ms' % (func.func_name, (t2-t1)*1000.0)
        return res
    return wrapper

# ====================================================================

def now_date_n_time() :
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
  return ss

# ====================================================================

# declare the @ decorator just before the function, invokes print_timing()
#@print_timing
def echo_and_execute_function(cmd, gs, ls, id = None) :
  print 'FUNCTION: |%s|' % cmd
  start_time = time.clock()
  if ( id is not None ) :
    pass
    print '   ID: %s' % id
  pass
  o = compile(cmd, '<string>', 'eval')
  rc = eval(o, gs, ls)
  stop_time = time.clock()
  delta_time = stop_time - start_time
  print '%s) rc =|%s|   DELTA TIME = %f sec (NOW: %s)' % (id, rc, delta_time, now_date_n_time())
  return rc

# ====================================================================

# declare the @ decorator just before the function, invokes print_timing()
@print_timing
def echo_and_execute_cmd(cmd, id = None) :
  print 'COMMAND: |%s|' % cmd
  start_time = time.clock()
  if ( id is not None ) :
    pass
    print '   ID: |%s|' % id
  pass
  rc = os.popen(cmd).read()
  stop_time = time.clock()
  delta_time = stop_time - start_time
  print '%s) rc =|%s|   DELTA TIME = %f sec (NOW: %s)' % (id, rc, delta_time, now_date_n_time())
  return rc

# ====================================================================

def start_function(id, cmd, gs, ls, mutexes, lock = None) :
  if( lock is not None ) : lock.acquire()
  mutexes[id] = id
  if( lock is not None ) : lock.release()
  rc = thread.start_new(start_function_thread, (id, cmd, gs, ls, mutexes, lock))
  return rc


# ====================================================================

def start_function_thread(id, cmd, gs, ls, mutexes, lock = None) :

  #lock.acquire()
  #mutexes[id] = id
  #lock.release()

  #print 'MUTEXES:' , mutexes
  #print 'BEFORE FUNCTION EXECUTION id = %s -=>' % id, mutexes
  start_time = time.clock()
  rc = echo_and_execute_function(cmd, gs, ls, id)
  stop_time = time.clock()
  delta_time = stop_time - start_time

  if( lock is not None ) : lock.acquire()
  mutexes[id] = 0  # signal waiter that this process is done
  if( lock is not None ) : lock.release()

  #print '> AFTER FUNCTION THREAD EXECUTION id = %s DELTA TIME = %f sec (NOW: %s)-=>' % (id, delta_time, now_date_n_time()), mutexes
  #print 'MUTEXES:' , mutexes
  #print 'THREAD |%s| IS COMPLETE' % id
  return rc

# ====================================================================

def start_command(id, cmd, mutexes, lock = None) :
  if( lock is not None ) : lock.acquire()
  mutexes[id] = id
  if( lock is not None ) : lock.release()
  rc = thread.start_new(start_command_thread, (id, cmd, mutexes, lock))
  return rc
  
# ====================================================================

def start_command_thread(id, cmd, mutexes, lock = None) :

  #lock.acquire()
  #mutexes[id] = id
  #lock.release()

  #print 'MUTEXES:' , mutexes
  #print 'BEFORE COMMAND EXECUTION id = %s -=>' % id, mutexes 
  start_time = time.clock()
  rc = echo_and_execute_cmd(cmd, id)
  stop_time = time.clock()
  delta_time = stop_time - start_time

  if( lock is not None ) : lock.acquire()
  mutexes[id] = 0  # signal waiter that this process is done
  if( lock is not None ) : lock.release()
  
  #print '* AFTER COMMAND THREAD EXECUTION id = %s DELTA TIME = %f sec (NOW: %s)-=>' % (id, delta_time, now_date_n_time()), mutexes
  #print 'MUTEXES:' , mutexes
  #print 'THREAD |%s| IS COMPLETE' % id
  return rc

# ====================================================================

def wait_for_mutexes_integer_ids(mutexes, cmd = '') :
  #print '*** WAITING FOR THREADS TO COMPLETE ***'
  last_sum = -1
  current_sum = sum(mutexes.values())
  while current_sum :
    last_sum = current_sum
    current_sum = sum(mutexes.values())
    if( last_sum != current_sum ) :
      pass
      #print 'THREAD ENDED: ', mutexes
    pass
  if ( not len(cmd) == 0 ) : 
    rc = echo_and_execute_cmd(cmd)
    #print 'FINAL: rc = %s' % rc
  pass
  # mutexes = {}

# ====================================================================

def wait_for_mutexes(mutexes, cmd = '') :
  #print '*** WAITING FOR THREADS TO COMPLETE ***'
  last_proc_list = list()
  current_proc_list = sorted(mutexes.values())
  while any(mutexes.values()) :
    last_proc_list = current_proc_list
    current_proc_list =  sorted(mutexes.values())
    # print '\n\n any(mutexes.values()) = ', any(mutexes.values())
    # print ' current_proc_list = ', current_proc_list
    # print ' last_proc_list = ', last_proc_list
    if( last_proc_list != current_proc_list ) :
      #su =  list(set(last_proc_list) - set(current_proc_list))
      su =  list((set(last_proc_list)).difference(set(current_proc_list)))
      #print 'THREAD ENDED: ', su
      #print 'MUTEXES NOW =->', mutexes
      su =  list(set(last_proc_list) - set(current_proc_list))
      #print '\n\n any(mutexes.values()) = ', any(mutexes.values())
      #print ' current_proc_list = ', current_proc_list
      #print ' last_proc_list = ', last_proc_list
    pass
  if ( not len(cmd) == 0 ) : 
    rc = echo_and_execute_cmd(cmd, 'FINAL')
    #print 'FINAL: rc = %s' % rc
  #print '\n\n any(mutexes.values()) = ', any(mutexes.values())
  #print ' current_proc_list = ', current_proc_list
  #print ' last_proc_list = ', last_proc_list
  #print ' NOW proc_list = ', sorted(mutexes.values())
  # mutexes = {}

# ====================================================================
# ID's DO NOT NEED TO BE INTEGERS ANY MORE BUT THEY STILL NEED TO BE UNIQUE

if __name__ == "__main__":
  
  def test_func(mess) :
    print 'TEST_FUNCTION: %s' % mess
    return '** ' + mess 

  thread_lock = thread.allocate_lock()
  #thread_lock = None
  mutexes = {}
  start_command('X Y', 'START', mutexes, thread_lock)
  start_function('XX', 'test_func(\'MESSAGE ONE\')', globals(), locals(), mutexes, thread_lock)

  start_command(3, 'START', mutexes, thread_lock)
  start_function(4, 'test_func(\'MESSAGE TWO\')', globals(), locals(), mutexes, thread_lock)

  start_command(5, 'START', mutexes, thread_lock)
  start_function(6, 'test_func(\'MESSAGE THREE\')', globals(), locals(), mutexes, thread_lock)

  start_command(7, 'START', mutexes)
  start_function(8, 'test_func(\'MESSAGE FOUR - 4\')', globals(), locals(), mutexes)
  
  wait_for_mutexes(mutexes, 'echo "ALL THREADS COMPLETED"')

  echo_and_execute_cmd('echo "ALL DONE"', '-=> DONE <=-')
