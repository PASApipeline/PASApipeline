import sys, os, subprocess, time




def run_commands (cmd_list, max_simult_processes=10):

    procs = []
    failed_procs = []
    
    for cmd in cmd_list:

        while len(procs) >= max_simult_processes:

            time.sleep(10)
            
            # wait until one opens up
            not_done = []

            for p in procs:
                ret = p['proc'].poll()
                if ret != None:
                    if ret != 0:
                        # errored out
                        failed_procs.append(p)
                        print >> sys.stderr, "process errored: CMD: " + p['cmd']
                    else:
                        # exited normally
                        print >> sys.stderr, "process completed successfully"
                else:
                    # process still running
                    not_done.append(p)
                    print >> sys.stderr, "process still running"

            procs = not_done
        
            
        # execute command
        print >> sys.stderr, "CMD: " + cmd

        p = { 'proc':None, 'cmd':None }
        p['proc'] = subprocess.Popen(cmd, shell=True)
        p['cmd'] = cmd
        procs.append(p)


    # wait for remaining ones to complete
    for p in procs:
        ret = p['proc'].wait()
       
        if ret != 0:
            # errored out
            failed_procs.append(p)
            print >> sys.stderr, "process errored: CMD: " + p['cmd'] 
        else:
            # exited normally
            print >> sys.stderr, "process completed successfully"



    print >> sys.stderr, "There were " + str(len(failed_procs)) + " failed processes."
    
    return failed_procs


if __name__ == "__main__":

    test_cmds = [ 'sleep 1',
                  'sleep 2',
                  'sleep 3',
                  'sleep 4',
                  'sleep x',
                  'sleep 6',
                  'doh',
                  'sleep 2' ]
                      
    run_commands(test_cmds, 2)
    
    sys.exit(0)
