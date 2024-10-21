###############################################################################
##                                                                           ##
## This code was written and is being maintained by Matias Villagra,         ##
## PhD Student in Operations Research @ Columbia, supervised by              ##
## Daniel Bienstock.                                                         ##
##                                                                           ##           
## For code readability, we make references to equations in:                 ##
## [1] D. Bienstock, and M. Villagra, Accurate Linear Cutting-Plane          ##
## Relaxations for ACOPF, arXiv:2312.04251v2, 2024                           ##
##                                                                           ##
## Please report any bugs or issues to: mjv2153@columbia.edu                 ## 
##                                                                           ##
## Jul 2024                                                                  ##
###############################################################################

import sys
import os
import numpy as np
import time
import reader
from myutils import *
from versioner import *
from log import danoLogger
from cutplane_mtp import gocutplane

def read_config(log, filename):

    log.joint("reading config file " + filename + "\n")

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.stateandquit("cannot open file", filename)
        sys.exit("failure")

    casefilename                 = 'NONE'
    lpfilename                   = ""
    lpfilename_cuts              = ""
    linenum                      = 0

    # Time-periods
    T                            = 4 # Only three options: 4, 12, 24
    
    # i2 definition
    i2                           = 0

    # Default rho parameter of i2(rho)+ c.f. Section 4.5 ineqs (19a), (19b)
    rho_threshold                = 1e2  
    
    # Convex non linear inequalities Jabr, i2 and limit
    jabr_inequalities            = 0
    i2_inequalities              = 0
    limit_inequalities           = 0

    # Cut management default parameters
    loud_cuts                    = 0 # Print cut information
    parallel_check               = 1 # Check whether new cut is parallel to incumbent cuts
    tolerance                    = 1e-5 # Target inequality violation threshold
    threshold_dotprod            = 5e-06 # Parallel check tolerance
    cut_age_limit                = 5
    tolerance                    = 1e-5 # Target violation threshold    
    jabrcuts                     = 1    # Jabr cuts active/inactive
    jabr_validity                = 0    # Check for cut validity
    dropjabr                     = 1    # Jabr cut management active/inactive
    most_violated_fraction_jabr  = 0.55 
    threshold                    = 1e-5 # Jabr ineq violation threshold        
    i2cuts                       = 0    # i2 cuts active/inactive
    i2_validity                  = 0    
    dropi2                       = 0
    most_violated_fraction_i2    = 0.15
    threshold_i2                 = 1e-3 # i2 ineq violation threshold    
    limitcuts                    = 0    # Limit cuts active/inactive
    limit_validity               = 0
    droplimit                    = 0
    most_violated_fraction_limit = 1
    threshold_limit              = 1e-5 # Limit ineq violation threshold
    
    # Solver default parameters
    solver_method                = 2
    crossover                    = 0
    barconvtol                   = 1e-6
    feastol                      = 1e-6
    opttol                       = 1e-6
    
    # Some cutting-plane algorithm parameters
    ftol                         = 1e-3
    ftol_iterates                = 5    
    max_time                     = 1200
    max_rounds                   = 100    

    # Previously computed solution and model output
    # default parameters
    fixflows                     = 0 # Fixes flows from loaded AC solution
    fixcs                        = 0 # Fixes c and s variable from loaded AC solution
    fix_tolerance                = 1e-5 # Tolerance for c,s and flow fixings 
    ampl_sol                     = 0 # Load AC solution
    writesol                     = 0 # Write cutting-plane solution
    writecuts                    = 0 # Write cuts to a .txt file
    addcuts                      = 0 # Load previously computed cuts
    writelps                     = 0 # Write cutting-plane linearly constrained models to .lp file
    writelastLP                  = 0 # Write the last iteration of the linearly constrained model to a .lp file

    # Dual information
    getduals                     = 0

    
    while linenum < len(lines):
        thisline = lines[linenum].split()
        if len(thisline) > 0:

            if thisline[0] == 'casefilename':
                casefilename = thisline[1]

            elif thisline[0] == 'lpfilename':
                lpfilename = thisline[1]

            elif thisline[0] == 'lpfilename_cuts':
                lpfilename_cuts = thisline[1]
            
            elif thisline[0] == 'jabr_inequalities':
                jabr_inequalities = 1

            elif thisline[0] == 'dropjabr':
                dropjabr     = 1

            elif thisline[0] == 'dropi2':
                dropi2 = 1                

            elif thisline[0] == 'limit_inequalities':
                limit_inequalities = 1

            elif thisline[0] == 'solver_method':
                solver_method = int(thisline[1])

            elif thisline[0] == 'i2cuts':
                i2cuts           = 1

            elif thisline[0] == 'most_violated_fraction_i2':
                most_violated_fraction_i2 = float(thisline[1])                
                
            elif thisline[0] == 'cut_age_limit':
                cut_age_limit = int(thisline[1])

            elif thisline[0] == 'threshold':
                threshold = float(thisline[1])

            elif thisline[0] == 'threshold_i2':
                threshold_i2 = float(thisline[1])

            elif thisline[0] == 'threshold_dotprod':
                threshold_dotprod = float(thisline[1])

            elif thisline[0] == 'parallel_check':
                parallel_check    = 1
                
            elif thisline[0] == 'tolerance':
                tolerance = float(thisline[1])

            elif thisline[0] == 'most_violated_fraction_jabr':
                most_violated_fraction_jabr = float(thisline[1])

            elif thisline[0] == 'crossover':
                crossover = int(thisline[1])

            elif thisline[0] == 'max_time':
                max_time  = float(thisline[1])

            elif thisline[0] == 'max_rounds':
                max_rounds = int(thisline[1])

            elif thisline[0] == 'i2_inequalities':
                i2_inequalities = 1

            elif thisline[0] == 'jabr_validity':
                jabr_validity = 1

            elif thisline[0] == 'i2_validity':
                i2_validity = 1

            elif thisline[0] == 'fixflows':
                fixflows = 1

            elif thisline[0] == 'ampl_sol':
                ampl_sol = 1

            elif thisline[0] == 'fixcs':
                fixcs    = 1

            elif thisline[0] == 'fix_tolerance':
                fix_tolerance = float(thisline[1])

            elif thisline[0] == 'most_violated_fraction_loss':
                most_violated_fraction_loss = float(thisline[1])

            elif thisline[0] == 'FeasibilityTol':
                FeasibilityTol = float(thisline[1])

            elif thisline[0] == 'writecuts':
                writecuts  = 1

            elif thisline[0] == 'addcuts':
                addcuts     = 1
                
            elif thisline[0] == 'droplimit':
                droplimit  = 1

            elif thisline[0] == 'most_violated_fraction_limit':
                most_violated_fraction_limit    = 1

            elif thisline[0] == 'threshold_limit':
                threshold_limit  = float(threshold_limit)

            elif thisline[0] == 'limit_validity':
                limit_validity  = 1

            elif thisline[0] == 'limitcuts':
                limitcuts    = 1

            elif thisline[0] == 'writelps':
                writelps     = 1

            elif thisline[0] == 'ftol':
                ftol          = float(thisline[1])

            elif thisline[0] == 'ftol_iterates':
                ftol_iterates = int(thisline[1])

            elif thisline[0] == 'jabrcuts':
                jabrcuts      = 1

            elif thisline[0] == 'loud_cuts':
                loud_cuts     = 1

            elif thisline[0] == 'i2':
                i2            = 1

            elif thisline[0] == 'writesol':
                writesol      = 1

            elif thisline[0] == 'T':
                T           = int(thisline[1])

            elif thisline[0] == 'barconvtol':
                barconvtol    = float(thisline[1])

            elif thisline[0] == 'feastol':
                feastol       = float(thisline[1])

            elif thisline[0] == 'opttol':
                opttol        = float(thisline[1])

            elif thisline[0] == 'writelastLP':
                writelastLP   = 1

            elif thisline[0] == 'rho_threshold':
                rho_threshold = float(thisline[1])

            elif thisline[0] == 'getduals':
                getduals = 1
                
            elif thisline[0] == 'END':
                break
                
            else:
                sys.exit("main_mtp: illegal input " + thisline[0] + " bye")


        linenum += 1

    all_data                        = {}
    all_data['casefilename']        = casefilename
    casename = all_data['casename'] = casefilename.split('data/')[1].split('.m')[0] 
    all_data['lpfilename']          = all_data['casename'] + '.lp'
    all_data['lpfilename_cuts']     = all_data['casename'] + '_cuts.lp'    

    if len(lpfilename):
        all_data['lpfilename']  = lpfilename
    
    if len(lpfilename_cuts):
        all_data['lpfilename_cuts'] = lpfilename_cuts

        
    # Convex non linear inequalities Jabr, i2 and limit
    all_data['jabr_inequalities']             = jabr_inequalities
    all_data['i2_inequalities']               = i2_inequalities
    all_data['limit_inequalities']            = limit_inequalities

    # Some cutplane parameters
    all_data['max_rounds']                    = max_rounds
    all_data['max_time']                      = max_time
    all_data['ftol']                          = ftol
    all_data['ftol_iterates']                 = ftol_iterates
    all_data['T']                             = T
    all_data['rho_threshold']                 = rho_threshold 

    # Previously computed solutions and options for writing
    # and loading cuts
    all_data['ampl_sol']                      = ampl_sol        
    all_data['fix_tolerance']                 = fix_tolerance
    all_data['fixflows']                      = fixflows
    all_data['fixcs']                         = fixcs
    all_data['writelastLP']                   = writelastLP    
    all_data['writelps']                      = writelps
    all_data['writesol']                      = writesol    
    all_data['writecuts']                     = writecuts
    all_data['addcuts']                       = addcuts

    # Solver parameters
    all_data['solver_method']                 = solver_method
    all_data['crossover']                     = crossover    
    all_data['barconvtol']                    = barconvtol
    all_data['feastol']                       = feastol
    all_data['opttol']                        = opttol

    # Dual information
    all_data['getduals']                      = getduals  
    all_data['duals']                         = {}
    all_data['dual_diff']                     = {}

    # Cuts
    all_data['threshold']              = threshold 
    all_data['tolerance']              = tolerance 
    all_data['threshold_dotprod']      = threshold_dotprod 
    all_data['cut_age_limit']          = cut_age_limit
    all_data['loud_cuts']              = loud_cuts
    all_data['parallel_check']         = parallel_check
    
    all_data['jabrcuts'] = jabrcuts
    if jabrcuts:
        all_data['most_violated_fraction_jabr'] = most_violated_fraction_jabr
        all_data['jabr_cuts']                   = {}
        all_data['num_jabr_cuts_rnd']           = {}
        all_data['num_jabr_cuts_added']         = 0
        all_data['num_jabr_cuts_dropped']       = 0
        all_data['dropped_jabrs']               = []    
        all_data['jabr_cuts_info']              = {}
        all_data['max_error_jabr']              = 0
        all_data['total_jabr_dropped']          = 0

    all_data['ID_jabr_cuts']         = 0 
    all_data['num_jabr_cuts']        = 0 
    all_data['dropjabr']             = dropjabr
    all_data['jabr_validity']        = jabr_validity

    all_data['limitcuts'] = limitcuts
    if limitcuts:
        all_data['most_violated_fraction_limit'] = most_violated_fraction_limit
        all_data['limit_cuts']                   = {}
        all_data['num_limit_cuts_rnd']           = {}
        all_data['num_limit_cuts_added']         = 0
        all_data['num_limit_cuts_dropped']       = 0
        all_data['limit_cuts_info']              = {}
        all_data['threshold_limit']              = threshold_limit
        all_data['dropped_limit']                = []
        all_data['max_error_limit']              = 0
        all_data['total_limit_dropped']          = 0

    all_data['ID_limit_cuts']   = 0 
    all_data['num_limit_cuts']  = 0 
    all_data['droplimit']       = droplimit
    all_data['limit_validity']  = limit_validity    


    all_data['i2cuts']                              = i2cuts
    all_data['i2']                                  = i2

    if (all_data['i2cuts'] == 1) or (all_data['i2_inequalities']):
        all_data['i2'] = 1
        
    if i2cuts:
        all_data['most_violated_fraction_i2'] = most_violated_fraction_i2
        all_data['i2_cuts']                   = {}
        all_data['num_i2_cuts_rnd']           = {}
        all_data['num_i2_cuts_added']         = 0
        all_data['num_i2_cuts_dropped']       = 0
        all_data['i2_cuts_info']              = {}
        all_data['threshold_i2']              = threshold_i2
        all_data['dropped_i2']                = []
        all_data['max_error_i2']              = 0
        all_data['total_i2_dropped']          = 0

    all_data['ID_i2_cuts']        = 0 
    all_data['num_i2_cuts']       = 0 
    all_data['dropi2']            = dropi2
    all_data['i2_validity']       = i2_validity

    ###########
    
    casetype = ''

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
        
    return all_data
        

if __name__ == '__main__':
    if len(sys.argv) > 4:
        print ('Usage: main_mtp.py file.conf [sols] [logfile]\n')
        exit(0)

    T0 = time.time()
    
    if len(sys.argv) == 4:
        sols      = sys.argv[2] + '/' # Optional subdirectory to where CP solutions can be written
        mylogfile = sys.argv[3]       # Optinal logfile
    else:
        sols      = ""
        mylogfile = "main.log"

    
    log = danoLogger(mylogfile)
        
    log.joint('\n  *********************************************************\n')
    log.joint(' ***********************************************************\n')
    log.joint(' ****                                                   ****\n')
    log.joint(' ****    Initializing AC-OPF Cutting-plane algorithm    ****\n')
    log.joint(' ****                                                   ****\n')
    log.joint(' ***********************************************************\n')
    log.joint('  *********************************************************\n\n')

    stateversion(log)
    
    all_data              = read_config(log,sys.argv[1])
    all_data['T0']        = T0
    all_data['sols']      = sols
    all_data['mylogfile'] = mylogfile
    all_data['datetime']  = mylogfile.strip("CPexp_").strip(".log")

    readcode = reader.readcase(log,all_data,all_data['casefilename'])

    gocutplane(log,all_data)
    
    log.closelog()
    
