################################################################################                            
##                                                                            ##                            
## This code was adapted and extended from the Gurobi OptiMod                 ##                            
## 'Optimal Power Flow', written by Daniel Bienstock. The extension was       ##                            
## written and it is being maintained by Matias Villagra                      ##                            
##                                                                            ##                            
## Please report any bugs or issues to                                        ##                            
##                        mjv2153@columbia.edu                                ##                            
##                                                                            ##                            
## Sept 2025                                                                  ##                            
################################################################################

import sys
import os
import time
from myutils import break_exit
from log import Logger
import math
from reader import *
from dc import *

def start_opf(logfilename):
    log = Logger(logfilename)

    alldata        = {}
    alldata['log'] = log
    alldata['LP']  = {}

    return alldata

    
def read_configfile(alldata, filename):
    """Function to read configurations for OPF solve from config file"""

    log = alldata['log']

    try:
        log.joint("Reading configuration file %s\n"%filename)
        f     = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.raise_exception("Error: Cannot open file %s\n"%filename)
        
    casefilename                 = 'NONE'
    dictionary_input             = False
    lpfilename                   = 'dcopf.lp'
    linenum                      = 0
    casetype                     = ''    
    writeLP                      = 0
    condnumber                   = 0

    # Read lines of configuration file and save options
    while linenum < len(lines):
        thisline = lines[linenum].split()

        if len(thisline) <= 0:
            linenum += 1
            continue

        if thisline[0][0] == '#':
            linenum += 1
            continue

        if thisline[0] == 'casefilename':
            casefilename = thisline[1]

        elif thisline[0] == 'dictionary_input':
            dictionary_input = True

        elif thisline[0] == 'lpfilename':
            lpfilename = thisline[1]

        elif thisline[0] == 'writeLP':
            writeLP = 1
            
        elif thisline[0] == 'condnumber':
            condnumber  = 1

        elif thisline[0] == 'END':
            break

        else:
            log.raise_exception("Error: Illegal input %s\n"%thisline[0])

        linenum += 1

    casename     = casefilename.split('data/')[1].split('.m')[0]
        
    log.joint("Settings:\n")
    for x in [('casefilename', casefilename), ('lpfilename', lpfilename),
              ('dictionary_input', dictionary_input),
              ('writeLP', writeLP),('casename',casename),
              ('condnumber',condnumber),('casetype',casetype)]:
        alldata[x[0]] = x[1]
        log.joint("  {} {}\n".format(x[0], x[1]))

    if alldata['casefilename'] == 'NONE':
       log.raise_exception('Error: No casefile provided\n')
       

if __name__ == '__main__':

    if len(sys.argv) > 3:
        log.raise_exception('Usage: main.py dc.conf [logfile]\n')

    T0 = time.time()
    
    if len(sys.argv) == 3:
        mylogfile = sys.argv[2] # Optinal logfile
    else:
        sols      = ""
        mylogfile = "dc.log"

    log = Logger(mylogfile)
    
    alldata = start_opf(mylogfile)
    log = alldata['log']

    alldata['t0'] = time.time()
    
    read_configfile(alldata, sys.argv[1])

    if alldata['dictionary_input']:
        main_dict = read_case_build_dict(alldata)
    else:
        read_case(alldata)
        

    lpformulator_dc(alldata)

    log.joint("\n writing casename, opt stauts, obj, " +
              "dinfs, runtime, numvars and constrs, and (kappa) to summary_dcopf.log\n")

    alldata['runtime'] = time.time() - alldata['t0']
    summary_dcopf   = open("summary_dcopf.log","a+") 

    if alldata['condnumber']:
        summary_dcopf.write(' case ' + alldata['casename'] + ' casetype '
                            + alldata['casetype'] + ' opt_status '
                            + str(alldata['optstatus']) + ' obj '
                            + str(alldata['objval']) + ' runtime '
                            + str(alldata['runtime']) + ' dinfs_sc '
                            + str(alldata['dinfs_scaled']) + ' dinfs '
                            + str(alldata['dinfs']) + ' numvars '
                            + str(alldata['numvars']) + ' numconstrs '
                            + str(alldata['numconstrs']) + ' kappa '
                            + str(alldata['Kappa'])
                            +  '\n')
        log.joint(' kappa ' + str(alldata['Kappa']) + ' kappaExact ' + str(alldata['KappaExact']) + '\n')

    else:
        summary_dcopf.write(' case ' + alldata['casename'] + ' casetype '
                            + alldata['casetype'] + ' opt_status '
                            + str(alldata['optstatus']) + ' obj '
                            + str(alldata['objval']) + ' runtime '
                            + str(alldata['runtime']) + ' dinfs_sc '
                            + str(alldata['dinfs_scaled']) + ' dinfs '
                            + str(alldata['dinfs']) + ' numvars '
                            + str(alldata['numvars']) + ' numconstrs '
                            + str(alldata['numconstrs']) + '\n')        
        

    summary_dcopf.close()    

    log.close_log()    
       
