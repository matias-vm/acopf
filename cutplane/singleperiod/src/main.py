###############################################################################
#
# This code was written and is being maintained by Matias Villagra,
# PhD Student in Operations Research @ Columbia, supervised by 
# Daniel Bienstock.
#
# Please report any bugs or issues (for sure there will be) to
#                         mjv2153@columbia.edu
#
# Oct 2023
###############################################################################

import sys
import os
import numpy as np
import time
import reader
from myutils import breakexit
from versioner import *
from log import danoLogger
from cutplane import gocutplane

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

    #i2 definition
    i2                           = 0
    
    # Convex non linear inequalities Jabr, i2 and limit 
    jabr_inequalities            = 0
    i2_inequalities              = 0
    limit_inequalities           = 0

    # Cut management default parameters
    writecuts                    = 0
    addcuts                      = 0    
    loud_cuts                    = 0
    tolerance                    = 1e-5
    jabrcuts                     = 1
    most_violated_fraction_jabr  = 0.55
    threshold                    = 1e-5
    dropjabr                     = 1
    jabr_validity                = 0
    i2cuts                       = 0    
    most_violated_fraction_i2    = 0.15
    dropi2                       = 0
    threshold_i2                 = 1e-03
    i2_validity                  = 0
    limitcuts                    = 0
    droplimit                    = 0    
    most_violated_fraction_limit = 1
    threshold_limit              = 1e-5    
    limit_validity               = 0
    parallel_check               = 0    
    threshold_dotprod            = 5e-06
    cut_age_limit                = 5

    # Solver default parameters
    solver_method                = 2
    crossover                    = 0

    # Some cutting-plane algorithm parameters                                                      
    ftol                         = 1e-5
    ftol_iterates                = 5
    max_time                     = 1000
    max_rounds                   = 100


    # Previously computed solution and model output                                                
    # default parameters
    ampl_sol                     = 0
    fixflows                     = 0
    fixcs                        = 0
    tol_fix                      = 1e-5
    writesol                     = 0
    writelps                     = 0    

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

            elif thisline[0] == 'fixcs':
                fixcs    = 1

            elif thisline[0] == 'tol_fix':
                tol_fix    = float(thisline[1])

            elif thisline[0] == 'most_violated_fraction_loss':
                most_violated_fraction_loss = float(thisline[1])

            elif thisline[0] == 'writecuts':
                writecuts  = 1

            elif thisline[0] == 'addcuts':
                addcuts    = 1

            elif thisline[0] == 'droplimit':
                droplimit  = 1

            elif thisline[0] == 'most_violated_fraction_limit':
                most_violated_fraction_limit    = 1

            elif thisline[0] == 'threshold_limit':
                threshold_limit     = float(threshold_limit)

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

            elif thisline[0] == 'ampl_sol':
                ampl_sol            = 1                

            elif thisline[0] == 'writesol':
                writesol      = 1

            elif thisline[0] == 'END':
                break
                
            else:
                sys.exit("illegal input " + thisline[0] + " bye")


        linenum += 1

    all_data                    = {}
    all_data['casefilename']    = casefilename
    all_data['casename']        = casefilename.split('data/')[1].split('.m')[0] 
    all_data['lpfilename']      = all_data['casename'] + '.lp'
    all_data['lpfilename_cuts'] = all_data['casename'] + '_cuts.lp'    

    if len(lpfilename):
        all_data['lpfilename']  = lpfilename
    
    if len(lpfilename_cuts):
        all_data['lpfilename_cuts'] = lpfilename_cuts
    

    all_data['jabr_inequalities']           = jabr_inequalities
    all_data['i2_inequalities']             = i2_inequalities
    all_data['limit_inequalities']          = limit_inequalities

    all_data['max_rounds']                    = max_rounds
    all_data['solver_method']                 = solver_method
    all_data['crossover']                     = crossover
    all_data['max_time']                      = max_time
    all_data['ftol']                          = ftol
    all_data['ftol_iterates']                 = ftol_iterates    
    
    #cuts
    all_data['threshold']                     = threshold
    all_data['tolerance']                     = tolerance
    all_data['threshold_dotprod']             = threshold_dotprod    
    all_data['cut_age_limit']                 = cut_age_limit
    all_data['parallel_check']                = parallel_check
    all_data['addcuts']                       = addcuts    
    all_data['writecuts']                     = writecuts
    all_data['writelps']                      = writelps
    all_data['loud_cuts']                     = loud_cuts

    
    all_data['jabrcuts'] = jabrcuts
    if jabrcuts:
        all_data['most_violated_fraction_jabr'] = most_violated_fraction_jabr
        all_data['ID_jabr_cuts']                = 0
        all_data['jabr_cuts']                   = {}
        all_data['num_jabr_cuts']               = 0
        all_data['num_jabr_cuts_rnd']           = {}
        all_data['num_jabr_cuts_added']         = 0
        all_data['num_jabr_cuts_dropped']       = 0
        all_data['dropped_jabrs']               = []    
        all_data['jabr_cuts_info']              = {}
        all_data['max_error_jabr']              = 0
        all_data['total_jabr_dropped']          = 0

    all_data['dropjabr']               = dropjabr
    all_data['jabr_validity']          = jabr_validity

    all_data['limitcuts'] = limitcuts
    if limitcuts:
        all_data['most_violated_fraction_limit']    = most_violated_fraction_limit
        all_data['ID_limit_cuts']                   = 0
        all_data['limit_cuts']                      = {}
        all_data['num_limit_cuts']                  = 0
        all_data['num_limit_cuts_rnd']              = {}
        all_data['num_limit_cuts_added']            = 0
        all_data['num_limit_cuts_dropped']          = 0
        all_data['limit_cuts_info']                 = {}
        all_data['threshold_limit']                 = threshold_limit
        all_data['dropped_limit']                   = []
        all_data['max_error_limit']                 = 0
        all_data['total_limit_dropped']             = 0

    all_data['droplimit']                       = droplimit
    all_data['limit_validity']                  = limit_validity    


    all_data['i2cuts']                              = i2cuts
    all_data['i2']                                  = i2

    if (all_data['i2cuts'] == 1) or (all_data['i2_inequalities']):
        all_data['i2'] = 1
        

    if i2cuts:
        all_data['most_violated_fraction_i2'] = most_violated_fraction_i2
        all_data['ID_i2_cuts']                = 0
        all_data['i2_cuts']                   = {}
        all_data['num_i2_cuts']               = 0
        all_data['num_i2_cuts_rnd']           = {}
        all_data['num_i2_cuts_added']         = 0
        all_data['num_i2_cuts_dropped']       = 0
        all_data['i2_cuts_info']              = {}
        all_data['threshold_i2']              = threshold_i2
        all_data['dropped_i2']                = []
        all_data['max_error_i2']              = 0
        all_data['total_i2_dropped']          = 0

    all_data['dropi2']                    = dropi2
    all_data['i2_validity']               = i2_validity
        
    all_data['ampl_sol']                      = ampl_sol    
    all_data['fixflows']                      = fixflows
    all_data['fixcs']                         = fixcs
    all_data['tol_fix']                       = tol_fix
    all_data['writesol']                      = writesol

    return all_data
        

if __name__ == '__main__':
    if len(sys.argv) > 3 or len(sys.argv) < 2:
        print ('Usage: main.py file.config [logfile]\n')
        exit(0)

    T0 = time.time()
    mylogfile = "main.log"

    if len(sys.argv) == 3:
        mylogfile = sys.argv[2]

    log = danoLogger(mylogfile)

    log.joint('\n  *********************************************************\n')
    log.joint(' ***********************************************************\n')
    log.joint(' ****                                                   ****\n')
    log.joint(' ****    Initializing AC-OPF Cutting-plane algorithm    ****\n')
    log.joint(' ****                                                   ****\n')
    log.joint(' ***********************************************************\n')
    log.joint('  *********************************************************\n\n')


    stateversion(log)

    all_data       = read_config(log,sys.argv[1])    
    all_data['T0'] = T0

    readcode = reader.readcase(log,all_data,all_data['casefilename'])

    gocutplane(log,all_data)
    
    log.closelog()
    
