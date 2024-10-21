###############################################################################
##                                                                           ##
## This code was written and is being maintained by Matias Villagra,         ##
## PhD Student in Operations Research @ Columbia, supervised by              ##
## Daniel Bienstock.                                                         ##
##                                                                           ##
## Please report any bugs or issues (for sure there will be) to              ##
##                        mjv2153@columbia.edu                               ##
##                                                                           ##
## April 2024                                                                ##
###############################################################################

import sys
import os
import time
import reader
from myutils import breakexit
from versioner import *
from log import danoLogger
from ac_mtp2 import *
from ac_mtp import *
from socp_mtp import *

def read_config(log, filename):

    log.joint("reading config file " + filename + "\n")

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.stateandquit("cannot open file", filename)
        exit(1)

    if filename == "ac_mtp.conf":
        ac   = 1
        socp = 0
    elif filename == "socp_mtp.conf":
        socp = 1
        ac   = 0
    else:
        log.joint("Please check name of config file\n")
        exit(0)
        
    casefilename       = '../data/none'
    T                  = 4                  # Default number of time-periods
    modfile            = '../modfiles/none' 
    solver             = 'knitroampl' # Default solver
    max_time           = "8000"       # Default solver time limit

    # i2plus or Jabr relaxation
    jabr   = 0 
    i2plus = 0
    
    # Use a previously computed AC solution
    initial_solution   = 0 # An AC solution is used to initialized the nonlinear solver
    initial_solution_0 = 0 # The first single-period solution is propagated and used to initialized the nonlinear solver
    ampl_sol           = 0 # Loads a previously computed AC solution

    # Knitro parameters
    feastol_abs          = "1e-6" # Absolute feasibility tolerance
    feastol_rel          = "1e-6" # Relative feasibility tolerance
    opttol_abs           = "1e-6" # Absolute optimality tolerance
    opttol_rel           = "1e-6" # Relative optimality tolerance
    scale                = "0"    # Scaling 
    ftol                 = "1e-5" # Threshold for objective relative improvement 
    ftol_iters           = "3"    # Iterations without sufficient objective improvement
    honorbnds            = "1"    # Satisfy/honor variable bounds
    linsolver            = "6"    # Type of linsolver
    linsolver_numthreads = "10"   # Number of threads linsolver
    blas_numthreads      = "1"    # Number of threads to use for the MKL BLAS
    blasoption           = "1"    # Specifies the BLAS/LAPACK function library
    wstart               = "0"    # Invoke a warm-start strategy 
    bar_initmu           = "1e-1" # Specifies the initial value for the barrier parameter
    bar_murule           = "0"    # Indicates which strategy to use for modifying the barrier parameter mu
    multistart           = 0      # Multistart

    # Write solution to .txt and .sol files
    writesol             = 0

    # Write model a human-readable format (similar in spirit to .lp file)
    expand               = 0

    # AC heuristics
    heuristic1           = 0
    heuristic2           = 0

    
    linenum = 0
    while linenum < len(lines):
        thisline = lines[linenum].split()
        if len(thisline) > 0:

            if thisline[0] == 'casefilename':
                casefilename  = thisline[1]

            elif thisline[0] == 'modfile':
                modfile       = thisline[1]
                
            elif thisline[0] == 'heuristic1':
                heuristic1   = 1
                heuristic2   = 0

            elif thisline[0] == 'heuristic2':
                heuristic1   = 0
                heuristic2   = 1                                
            
            elif thisline[0] == 'solver':
                solver        = thisline[1]

            elif thisline[0] == 'initial_solution':
                initial_solution = 1

            elif thisline[0] == 'initial_solution_0':
                initial_solution_0 = 1                

            elif thisline[0] == 'multistart':
                multistart         = 1

            elif thisline[0] == 'writesol':
                writesol           = 1

            elif thisline[0] == 'T':
                T                  = int(thisline[1])

            elif thisline[0] == 'expand':
                expand           = 1

            elif thisline[0] == 'feastol_abs':
                feastol_abs      = thisline[1]

            elif thisline[0] == 'opttol_abs':
                opttol_abs       = thisline[1]

            elif thisline[0] == 'feastol_rel':
                feastol_rel      = thisline[1]

            elif thisline[0] == 'opttol_rel':
                opttol_rel       = thisline[1]                

            elif thisline[0] == 'linsolver':
                linsolver    = thisline[1]

            elif thisline[0] == 'max_time':
                max_time     = thisline[1]

            elif thisline[0] == 'wstart':
                wstart       = '1'

            elif thisline[0] == 'bar_initmu':
                bar_initmu   = thisline[1]

            elif thisline[0] == 'ampl_sol':
                ampl_sol   = thisline[1]

            elif thisline[0] == 'i2plus':
                i2plus   = 1
                jabr     = 0
                
            elif thisline[0] == 'jabr':
                jabr     = 1
                i2plus   = 0
                
            elif thisline[0] == 'END':
                break

            else:
                exit("main_mtp: illegal input " + thisline[0] + "bye")

        linenum += 1


    all_data                         = {}
    all_data['casefilename']         = casefilename
    all_data['casename']             = casefilename.split('data/')[1].split('.m')[0]
    all_data['modfile']              = modfile.split('../modfiles/')[1]
    all_data['solver']               = solver
    all_data['initial_solution']     = initial_solution
    all_data['initial_solution_0']   = initial_solution_0    
    all_data['bar_initmu']           = bar_initmu
    all_data['multistart']           = multistart
    all_data['writesol']             = writesol
    all_data['T']                    = T
    all_data['expand']               = expand
    all_data['feastol_abs']          = feastol_abs
    all_data['opttol_abs']           = opttol_abs
    all_data['feastol_rel']          = feastol_rel
    all_data['opttol_rel']           = opttol_rel
    all_data['scale']                = scale
    all_data['honorbnds']            = honorbnds
    all_data['linsolver_numthreads'] = linsolver_numthreads
    all_data['blasoption']           = blasoption
    all_data['blas_numthreads']      = blas_numthreads
    all_data['ftol']                 = ftol
    all_data['ftol_iters']           = ftol_iters
    all_data['bar_murule']           = bar_murule
    all_data['linsolver']            = linsolver
    all_data['max_time']             = max_time
    all_data['wstart']               = wstart
    all_data['ac']                   = ac
    all_data['socp']                 = socp
    all_data['heuristic1']           = heuristic1
    all_data['heuristic2']           = heuristic2    
    all_data['ampl_sol']             = ampl_sol
    all_data['i2plus']               = i2plus
    all_data['jabr']                 = jabr

    
    casetype = ''
    casename = all_data['casename']
    
    if T == 4:
        all_data['uniform_drift'] = 0.01
        casetype = 'uniform_' + str(all_data['uniform_drift'])
        log.joint('case ' + all_data['casename'] + ' multi-period ' + str(T)
                  + ' ' + casetype + '\n')
        loadsfilename = '../data/mtploads/' + casename + '_mtploads_' + str(T) + '_u_' + str(all_data['uniform_drift']) + '.txt'
    elif T == 12:
        all_data['uniform_drift'] = 0.02
        casetype = 'uniform6_' + str(all_data['uniform_drift'])
        log.joint('case ' + all_data['casename'] + ' multi-period ' + str(T)
                  + ' ' + casetype + '\n')
        loadsfilename = '../data/mtploads/' + casename + '_mtploads_' + str(T) + '_u6_' + str(all_data['uniform_drift']) + '.txt'
    elif T == 24:
        all_data['uniform_drift'] = 0.01
        casetype = 'uniform5' + str(all_data['uniform_drift'])
        log.joint('case ' + all_data['casename'] + ' multi-period ' + str(T)
                  + ' ' + casetype + '\n')
        loadsfilename = '../data/mtploads/' + casename + '_mtploads_' + str(T) + '_u5_' + str(all_data['uniform_drift']) + '.txt'
    else:
        log.joint('Check config file, multi-period data only for T = 4, 12, and 24\n')
        exit(1)

    all_data['casetype']      = casetype
    all_data['loadsfilename'] = loadsfilename
    all_data['rampfilename']  = '../data/ramprates/' + casename + '_rampr_' + str(T) + '.txt'


    # Trimming casename for table

    if 'pglib' in casename:
        csplit = casename.split('_')
        rmcase = csplit[2].split('case')[1]
        casename       = rmcase + csplit[3] + '-' + csplit[5]
        casename_table = casename.split('.')[0]
    elif 'ACTIVS' in casename:
        csplit         = casename.split('_')
        casename_table = csplit[1].split('.')[0]
    else:
        casename_table = casename.split('_')[0].split('case')[1].split('.')[0]

    all_data['casename_table'] = casename_table
    

    return all_data

if __name__ == '__main__':
    if len(sys.argv) > 4:
        print ('Usage: main_mtp.py file.conf [sols] [logfile]\n')
        exit(1)

    T0 = time.time()

    if len(sys.argv) == 4: 
        sols      = sys.argv[2] + '/' # Optional subdirectory to where solutions can be written
        mylogfile = sys.argv[3]       # Optinal logfile 
    else:
        sols      = ""        
        mylogfile = "main.log"

    
    log = danoLogger(mylogfile)
    stateversion(log)

    all_data              = read_config(log,sys.argv[1])
    all_data['T0']        = T0
    all_data['sols']      = sols
    all_data['mylogfile'] = mylogfile
    all_data['datetime']  = mylogfile.strip("exp_").strip(".log")
    
    readcode = reader.readcase(log,all_data,all_data['casefilename'])

    T        = all_data['T']
    casename = all_data['casename']
    casetype = all_data['casetype']
    datetime = all_data['datetime']

    
    if all_data['ac']:
        all_data['solver']     = "knitroampl"
        #all_data['opttol_rel'] = "1e-3"
        #all_data['opttol_abs'] = "1e20"

        all_data['h1'] = {}
        all_data['h1']['objvalue'] = -1
        all_data['h1']['time']     = 0
        all_data['h1']['status']   = '-'
        all_data['h2'] = {}
        all_data['h2']['objvalue'] = -1
        all_data['h2']['time']     = 0
        all_data['h2']['status']   = '-'        
        
        
        log.joint('\n==================================================='
                  + '===========================\n')
        log.joint('===================================================='
                  + '==========================\n')        
        log.joint('\n                  Initializing Multi-Time Period ACOPF\n')
        if all_data['heuristic1']:
            all_data['modfile'] = 'ac_mtp_definedvars.mod'
            log.joint('\n**Running Heuristic 1: We give to Knitro the '
                      + 'full Multi-Time Period Formulation\n')
            log.joint(' feastol (rel/abs) = ' + str(all_data['feastol_rel']) + "/"
                      + str(all_data['feastol_abs'])
                      + ' opttol (rel/abs) = ' + str(all_data['opttol_rel']) + "/"
                      + str(all_data['opttol_abs']) 
                      + '\n')                        
            retcode = goac_mtp(log,all_data)
        elif all_data['heuristic2']:
            all_data['modfile'] = 'ac_definedvars.mod'
            all_data['max_time']     = "1200" 
            all_data['feastol_rel']  = "1e-3" # Hardcoded feasibility tolerances
            all_data['feastol_abs']  = "1e20"            
            log.joint('\n\nRunning Heuristic 2: We give to Knitro '
                      + 'one period at a time,\n')
            log.joint('and impose the ramping constraints using '
                      + 'generation from t-1\n')
            log.joint(' feastol (rel/abs) = ' + str(all_data['feastol_rel']) + "/"
                      + str(all_data['feastol_abs'])
                      + ' opttol (rel/abs) = ' + str(all_data['opttol_rel']) + "/"
                      + str(all_data['opttol_abs']) + "\n") 
            log.joint (' max time (secs) per period ' + str(all_data['max_time'])
                      + '\n')                        
            retcode = goac_mtp2(log,all_data)            
               
        # Set as default option when running MTP ACOPF
        else:
            
            log.joint('\n**Running Heuristic 3: First we run H1. If '
                      + 'it fails, relax tolerances and run H2.\n')
            log.joint(' feastol (rel/abs) = ' + str(all_data['feastol_rel']) + "/"
                      + str(all_data['feastol_abs'])
                      + ' opttol (rel/abs) = ' + str(all_data['opttol_rel']) + "/"
                      + str(all_data['opttol_abs']) + "\n")
            log.joint("max time (secs) per period " + str(all_data['max_time'])
                      + "\n")
                      
            all_data['modfile'] = 'ac_mtp_definedvars.mod'
            retcode  = goac_mtp(log,all_data)
            if retcode != 0:
                
                log.joint('\n\nRunning Heuristic 2: We give to Knitro '
                          + 'one period at a time,\n')
                log.joint('and impose the ramping constraints using '
                          + 'generation from t-1\n')
                log.joint(' We relax feasibility and optimality tolerances\n')
                all_data['modfile']      = 'ac_definedvars.mod'
                all_data['max_time']     = "1800" #| for T=12 1800 | T=4 1200 | T=24 4000
                all_data['feastol_rel']  = "1e-3" # Hardcoded feasibility tolerances
                all_data['feastol_abs']  = "1e20"
                all_data['bar_murule']   = "1"
                log.joint(' feastol (rel/abs) = ' + str(all_data['feastol_rel']) + "/"
                          + str(all_data['feastol_abs'])
                          + ' opttol (rel/abs) = ' + str(all_data['opttol_rel']) + "/"
                          + str(all_data['opttol_abs']) 
                          + ' max time (secs) per period ' + str(all_data['max_time'])
                          + '\n')
                retcode = goac_mtp2(log,all_data)

            # Writing results to a table
            log.joint(' writing results to table_ac\n')
            tablename_ac  = 'table_ac_' + datetime + '.txt' 
            table_ac      = open(tablename_ac,"a")
            casename_table   = all_data['casename_table']
            numvars    = all_data['numvars']
            numconstrs = all_data['numconstrs']

            if all_data['h1']['status'] == "0":
                objval = all_data['h1']['objvalue']
                time   = all_data['h1']['time']
                result = casename_table + " & " + str(numvars) + " & " + str(numconstrs) + " & " + str(objval) + " & " + str(time) + " &\\\ \n"
                table_ac.write(result)
            else:
                objval = all_data['h2']['objvalue']
                time   = all_data['h2']['time']
                status = all_data['h2']['status']
                if all_data['h2']['status'] == "0":
                    result = casename_table + " & " + str(numvars) + " & " + str(numconstrs) + " & " + str(objval) + " & " + str(time) + " &\\\ \n"
                else:
                    result = casename_table + " & " + str(numvars) + " & " + str(numconstrs) + " & " + status + " & " + str(time) + " &\\\ \n"
                table_ac.write(result)
                
            table_ac.close()
                
    elif all_data['socp']:
        log.joint('\n==================================================='
                  + '===========================\n')
        log.joint('==================================================='
                  + '===========================\n')        
        log.joint('\n      Initializing the JABR SOCP relaxation for '
                  + 'Multi-Time Period ACOPF\n\n')

        all_data['honorbnds']  = "1"
        all_data['scale']      = "0"                    
        
        if all_data['i2plus']:
            all_data['modfile'] = 'i2plus_mosek_mtp.mod'
            retcode             = i2plus_mosek_mtp(log,all_data)

        elif all_data['jabr']:
            if all_data['solver'] == 'mosek':          
                all_data['modfile'] = 'jabr_mosek_mtp.mod'
                retcode             = gosocp_mosek_mtp(log,all_data)
            else:
                all_data['modfile'] = 'jabr_mtp.mod'
                retcode             = gosocp_mtp(log,all_data)
        else:
            log.joint('main_mtp: i2plus or Jabr relxation? Please check config file.\n')
            exit(1)            

    else:
        log.joint(' modfile chosen ' + all_data['modfile'] + '\n')
        log.joint(' solver chosen ' + all_data['solver'] + '\n')
        log.joint('main_mtp: ac or socp? Please check config file.\n')
        exit(1)
        
    log.closelog()

