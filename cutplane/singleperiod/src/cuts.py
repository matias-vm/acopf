###############################################################################
##                                                                           ##
## This code was written and is being maintained by Matias Villagra,         ##
## PhD Student in Operations Research @ Columbia, supervised by              ##
## Daniel Bienstock.                                                         ##
##                                                                           ##
## Please report any bugs or issues (for sure there will be) to              ##
##                         mjv2153@columbia.edu                              ##
##                                                                           ##
## Oct 2023                                                                  ##
###############################################################################

from myutils import breakexit
from log import danoLogger
import time
import math
from json import dumps
from gurobipy import *
import numpy as np
from numpy import linalg as LA


def compute_normal(*args):
    n = len(args)
    v = np.zeros(n,dtype = 'float')
    for i in range(n):
        v[i] = args[i]

    norm = LA.norm(v)

    if norm == 0:
        log.joint(' BUG: check vector (cut)\n')
        brekaexit('check cut')
    return v / norm


def i2_cuts(log,all_data):
        
    log.joint('\n')
    log.joint(' **** i2-cuts ****\n')

    themodel       = all_data['themodel']
    IDtoCountmap   = all_data['IDtoCountmap']
    buses          = all_data['buses']
    branches       = all_data['branches']
    FeasibilityTol = all_data['tolerance']
    
    cvar      = all_data['cvar']
    Pvar_f    = all_data['Pvar_f']
    Qvar_f    = all_data['Qvar_f']
    i2var_f   = all_data['i2var_f']
    
    Pfvalues  = all_data['Pfvalues']
    Qfvalues  = all_data['Qfvalues']
    cvalues   = all_data['cvalues']
    i2fvalues = all_data['i2fvalues']

    rnd                      = all_data['round']
    i2_cuts                  = all_data['i2_cuts']
    i2_cuts_info             = all_data['i2_cuts_info']
    num_cuts                 = all_data['ID_i2_cuts']
    threshold                = all_data['threshold_i2']
    violated                 = {}
    violated_count           = 0
    num_i2_cuts_added        = 0
        
    all_data['NO_i2_cuts_violated'] = 0
    
    log.joint(' checking for violations of i2 inequalities ... \n')
    for branch in branches.values():
        
        f = branch.f
        t = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        violation = Pfvalues[branch]*Pfvalues[branch] + Qfvalues[branch]*Qfvalues[branch] - cvalues[buses[count_of_f]] * i2fvalues[branch]
        #violation_f = Pfvalues[branch]*Pfvalues[branch] + Qfvalues[branch]*Qfvalues[branch] - cvalues[buses[count_of_f]] * i2fvalues[branch]
        #violation_t = Ptvalues[branch]*Ptvalues[branch] + Qtvalues[branch]*Qtvalues[branch] - cvalues[buses[count_of_t]] * i2tvalues[branch]
        #violation = max(violation_f,violation_t)
        if violation > threshold:
            if all_data['loud_cuts']:
                log.joint('  violation ' + str(violation) + ' at branch ' 
                          + str(branch.count) + ' f ' + str(f) + ' t ' 
                          + str(t) + ' i2 value ' + str(i2fvalues[branch]) 
                          + '\n')
            violated_count  += 1
            violated[branch] = violation

    if violated_count == 0:
        all_data['NO_i2_cuts_violated'] = 1
        log.joint(' all i2 violations below threshold\n' )
        log.joint(' no more i2 cuts to add for current threshold\n' )
        return None

    log.joint('  number violated i2 ineqs ' + str(violated_count) + '\n')

    log.joint(' sorting most violated i2-envelope cuts ...\n')

    num_selected        =  math.ceil(violated_count * all_data['most_violated_fraction_i2'] )
    most_violated       = dict(sorted(violated.items(), key = lambda x: x[1], reverse = True)[:num_selected])
    most_violated_count = 0

    log.joint(' computing i2-envelope cuts ... \n')

    for branch in most_violated.keys():
        if (most_violated_count == 0):
            most_violated_branch = branch.count
            max_error = most_violated[branch]
                
        f                    = branch.f
        t                    = branch.t
        count_of_f           = IDtoCountmap[f]
        count_of_t           = IDtoCountmap[t]
        Pft                  = Pfvalues[branch]
        Qft                  = Qfvalues[branch]
        cff                  = cvalues[buses[count_of_f]]
        i2ft                 = i2fvalues[branch]

        cutnorm    = math.sqrt( (2 * Pft)**2 + (2 * Qft)**2 + (cff - i2ft)**2 )
        coeff_Pft  = 4 * Pft
        coeff_Qft  = 4 * Qft
        coeff_cff  = cff - i2ft - cutnorm
        coeff_i2ft = - (cff - i2ft) - cutnorm

        if all_data['parallel_check']:
            if parallel_check_i2(log,all_data,branch,coeff_Pft,coeff_Qft,
                                 coeff_cff,coeff_i2ft):
                continue

        most_violated_count += 1
        violation            = most_violated[branch]
        cutid                = num_cuts + most_violated_count

        if all_data['loud_cuts']:
            log.joint('  --> new i2-cut\n')
            log.joint('  branch ' + str(branch.count) + ' f ' + str(f) + ' t ' 
                      + str(t) + ' violation ' + str(violation) + ' cut id ' 
                      + str((num_cuts + most_violated_count)) + '\n' )
            log.joint('  values ' + ' Pft ' + str(Pft) + ' Qft ' + str(Qft) 
                      + ' cff ' + str(cff) + ' i2ft ' + str(i2ft) + '\n' )
            log.joint('  LHS coeff ' + ' Pft ' + str(coeff_Pft) + ' Qft ' 
                      + str(coeff_Qft) + ' cff ' + str(coeff_cff) + ' i2ft ' 
                      + str(coeff_i2ft) + '\n' )
            log.joint('  cutnorm ' + str(cutnorm) + '\n')

        #sanity check
        if all_data['i2_validity']:
            
            sol_Pf        = all_data['sol_Pfvalues'][branch]
            sol_Qf        = all_data['sol_Qfvalues'][branch]
            sol_c         = all_data['sol_cvalues'][branch]
            sol_s         = all_data['sol_svalues'][branch]
            sol_cbusf     = all_data['sol_cvalues'][buses[count_of_f]]
            sol_cbust     = all_data['sol_cvalues'][buses[count_of_t]]
            sol_i2f       = computei2value(log,all_data,branch,sol_c,sol_s,sol_cbusf,sol_cbust)
            violation     = coeff_Pft * sol_Pf + coeff_Qft * sol_Qf + coeff_cff * sol_cbusf + coeff_i2ft * sol_i2f
            relviolation  = violation / ( ( coeff_Pft**2 + coeff_Qft**2 + coeff_cff**2 + coeff_i2ft**2 )**0.5 ) 

            if relviolation > FeasibilityTol:
                log.joint('  WARNING, the i2 inequality associated to branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' is violated by the AC solution!\n')
                log.joint('  violation ' + str(violation) + '\n')
                log.joint('  relative violation ' + str(relviolation) + '\n')
                log.joint('  values (AC solution) ' + ' Pft ' + str(sol_Pf) + ' Qft ' + str(sol_Qf) + ' cff ' + str(sol_cbusf) + ' i2ft ' + str(sol_i2f) + '\n' )
                breakexit('check!')
            else:
                log.joint('  AC solution satisfies i2 inequality at branch ' + str(branch.count) + ' with slack ' + str(violation) + '\n')
        
        
        i2_cuts[(cutid,branch.count)] = (rnd,threshold)
        i2_cuts_info[branch][cutid]   = (rnd,violation,coeff_Pft,coeff_Qft,coeff_cff,coeff_i2ft,threshold,cutid)

        
        cutexp     = LinExpr()
        constrname = "i2_cut_"+str(cutid)+"_"+str(branch.count)+"r_"+str(rnd)+"_"+str(f)+"_"+str(t)
        cutexp    += coeff_Pft * Pvar_f[branch] + coeff_Qft * Qvar_f[branch] + coeff_cff * cvar[buses[count_of_f]] + coeff_i2ft * i2var_f[branch]

        themodel.addConstr(cutexp <= 0, name = constrname)

    log.joint('  number i2-envelope cuts added ' + str(most_violated_count) + '\n')
    log.joint('  max error (abs) ' + str(max_error) + ' at branch ' + str(most_violated_branch) + '\n' )


    all_data['most_violated_branch_i2'] = most_violated_branch
    all_data['max_error_i2']            = max_error
    all_data['ID_i2_cuts']             += most_violated_count
    all_data['num_i2_cuts_added']       = most_violated_count
    all_data['num_i2_cuts']            += most_violated_count
    all_data['num_i2_cuts_rnd'][rnd]    = most_violated_count


    if all_data['dropi2']:
        if all_data['addcuts']:
            t0_drop = time.time()
            drop_i2(log,all_data)
            t1_drop = time.time()
            log.joint('  time spent on drop i2 ' + str(t1_drop - t0_drop) + '\n')
        elif all_data['round'] >= all_data['cut_age_limit']: ######
            t0_drop = time.time()
            drop_i2(log,all_data)
            t1_drop = time.time()
            log.joint('  time spent on drop i2 ' + str(t1_drop - t0_drop) + '\n')


def drop_i2(log,all_data):
    
    #log.joint('\n')
    #log.joint(' **** drop i2-envelope cuts ****\n')
    log.joint(' dropping old and slack i2-envelope cuts ...\n')

    themodel             = all_data['themodel']
    branches             = all_data['branches']
    buses                = all_data['buses']
    IDtoCountmap         = all_data['IDtoCountmap']
    current_rnd          = all_data['round']
    cut_age_limit        = all_data['cut_age_limit']
    i2_cuts              = all_data['i2_cuts']
    i2_cuts_info         = all_data['i2_cuts_info']
    dropped_i2           = all_data['dropped_i2']
    drop_i2              = []
    num_i2_cuts_dropped  = all_data['num_i2_cuts_dropped'] 

    for key in i2_cuts.keys():  
        cut     = i2_cuts[key] 
        cut_rnd = cut[0]
        cut_age = current_rnd - cut_rnd

        if cut_age <= cut_age_limit:
            continue 
        
        cutid         = key[0]
        branchid      = key[1]
        branch        = branches[branchid]
        f             = branch.f
        t             = branch.t
        count_of_f    = IDtoCountmap[f]
        count_of_t    = IDtoCountmap[t]
        cut_threshold = cut[1]

        constrname    = "i2_cut_"+str(cutid)+"_"+str(branchid)+"r_"+str(cut_rnd)+"_"+str(f)+"_"+str(t)
        constr        = themodel.getConstrByName(constrname)
        slack         = - constr.getAttr("slack")

        if ( slack < - cut_threshold ):
            drop_i2.append(key)
            themodel.remove(themodel.getConstrByName(constrname))
            if all_data['loud_cuts']:
                log.joint('  --> removed i2-cut\n')
                log.joint('  cut ' + str(key) + ' branch ' 
                          + str(branchid) + ' f ' + str(f) + ' t ' + str(t) 
                          + ' was removed from the model\n')
    
    num_drop_i2       = len(drop_i2)
    
    if num_drop_i2:
        all_data['num_i2_cuts_dropped'] = num_drop_i2
        all_data['num_i2_cuts']        -= num_drop_i2
        all_data['total_i2_dropped']   += num_drop_i2

        dropped_i2.extend(drop_i2)

        for key in drop_i2:
            i2_cuts.pop(key)
        log.joint('  the cuts in drop_i2 list were removed from dict i2_cuts\n' )

        ##### drop i2-envelope cuts from i2_cuts_info
        for key in drop_i2:
            cutid = key[0]
            branchid = key[1]
            cuts_branch = i2_cuts_info[branches[branchid]]
            cuts_branch.pop(cutid)
        log.joint('  cuts in drop_i2 list were removed from i2_cuts_info\n') 
    else:
        all_data['num_i2_cuts_dropped'] = 0
        log.joint('  no i2-envelope cuts were dropped this round\n')
    

def limit_cuts(log,all_data):
        
    log.joint('\n')
    log.joint(' **** Limit-cuts ****\n') 

    themodel         = all_data['themodel']
    IDtoCountmap     = all_data['IDtoCountmap']
    buses            = all_data['buses']
    branches         = all_data['branches']
    Pvar_f           = all_data['Pvar_f']
    Qvar_f           = all_data['Qvar_f']
    Pvar_t           = all_data['Pvar_t']
    Qvar_t           = all_data['Qvar_t']
    FeasibilityTol   = all_data['tolerance']
    
    Pfvalues = all_data['Pfvalues']
    Qfvalues = all_data['Qfvalues']
    Ptvalues = all_data['Ptvalues']
    Qtvalues = all_data['Qtvalues']
    
    rnd                   = all_data['round']
    limit_cuts            = all_data['limit_cuts']    
    limit_cuts_info       = all_data['limit_cuts_info']
    num_cuts              = all_data['ID_limit_cuts']
    num_limit_cuts_added  = 0
    threshold             = all_data['threshold']
    violated              = {}
    violated_count        = 0
    
    all_data['NO_limit_cuts_violated'] = 0
    
    log.joint(' checking for violations of limit inequalities ... \n')
    for branch in branches.values():
        violation_f = Pfvalues[branch]*Pfvalues[branch] + Qfvalues[branch]*Qfvalues[branch] - branch.limit**2
        violation_t = Ptvalues[branch]*Ptvalues[branch] + Qtvalues[branch]*Qtvalues[branch] - branch.limit**2
        f_and_t = 0
        if violation_f > threshold:
            f_and_t        += 1
            violated_count += 1
            violated[branch] = (violation_f,'f')
        elif violation_t > threshold:
            f_and_t        += 1
            violated_count += 1
            violated[branch] = (violation_t,'t')
        if f_and_t == 2:
            log.joint('  -from and to- limit inequalities violated at branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + '\n')  
            #breakexit('f and t violated')
        f_and_t = 0

    if violated_count == 0:
        all_data['NO_limit_cuts_violated'] = 1
        log.joint(' all limit violations below threshold\n' )
        log.joint(' no more limit-cuts to add for current threshold\n' )
        return None
    
    log.joint('  number violated limits ' + str(violated_count) + '\n')

    log.joint(' sorting most violated limit-envelope cuts ...\n')

    num_selected        =  math.ceil(violated_count * all_data['most_violated_fraction_limit'] )
    most_violated       = dict(sorted(violated.items(), key = lambda x: x[1][0], reverse = True)[:num_selected])
    most_violated_count = 0

    log.joint(' computing limit-envelope cuts ... \n')

    for branch in most_violated.keys():
        if (most_violated_count == 0):
            most_violated_branch = branch.count
            max_error = most_violated[branch][0]

        from_or_to = most_violated[branch][1]
    
        if from_or_to == 'f':
            Pval                  = Pfvalues[branch]
            Qval                  = Qfvalues[branch]
        elif from_or_to == 't':
            Pval                 = Ptvalues[branch]
            Qval                 = Qtvalues[branch]
        else:
            log.joint(' we have a bug\n')
            breakexit('look for bug')

        f          = branch.f
        t          = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]

        #inequality P^2 + Q^2 <= u^2; cutting-plane: coeff_P * P + coeff_Q * Q <= z
        #t0 := constant in (0,1) s.t. (t0 * Pval)^2 + (t0 * Qval)^2 = u^2
        #newPval = t0 * Pval, newQval = t0 * Qval; (normal) separating hyperplane at newPval,newQval : (2 * newPval, 2 * newQval) 
        #Hence, coeff_P = 2 * newPval, coeff_Q = 2 * newQval, z = 2 * ( newPval^2 + newQval^2) = 2 * t0^2 * ( Pval^2 + Qval^2 ) = 2 * u^2

        violation  = most_violated[branch][0]
        u          = branch.limit
        u2         = u**2

        if violation + u2 < 1e-05:
            log.joint('  check branch\n') #case were U ~ 0 and P,Q > 0 (small)  if viol > threshold, this shouldnt happen (given threshold >= 1e-5)
            #breakexit('check')  

        #t0_sqred  = u**2 / ( violation + u**2 )
        #t0        = math.sqrt(t0_sqred)
        
        #t0         = u / math.sqrt(violation + u2)
        #coeff_P    = t0 * Pval # or coeff_P' = coeff_P / z and use 1 as RHS
        #coeff_Q    = t0 * Qval
        #z          = u2

        t0         = 1 / (u * math.sqrt(violation + u2))
        coeff_P    = t0 * Pval #RHS is 1
        coeff_Q    = t0 * Qval
        z          = 1

        if all_data['parallel_check']:
            if parallel_check_limit(log,all_data,branch,coeff_P,coeff_Q,
                                    from_or_to):
                continue

        #we add the cut
        most_violated_count += 1
        cutid                = num_cuts + most_violated_count

        if all_data['loud_cuts']:
            log.joint('  --> new cut\n')
            log.joint('  branch ' + str(branch.count) + ' f ' + str(f) + ' t ' 
                      + str(t) + ' violation ' + str(violation) + ' cut id ' 
                      + str(cutid) + '\n')
            if from_or_to == 'f':
                log.joint('  values ' + ' Pft ' + str(Pval) + ' Qft ' 
                          + str(Qval) + '\n')
                log.joint('  LHS coeff ' + ' Pft ' + str(coeff_P) + ' Qft ' 
                          + str(coeff_Q) + ' RHS ' + str(z) + '\n')
            elif from_or_to == 't':
                log.joint('  values ' + ' Ptf ' + str(Pval) + ' Qtf ' 
                          + str(Qval) + '\n')
                log.joint('  LHS coeff ' + ' Ptf ' + str(coeff_P) + ' Qtf ' 
                          + str(coeff_Q) + ' RHS ' + str(z) + '\n')
        
        #sanity check
        if all_data['limit_validity']:
            
            if from_or_to == 'f':
                sol_Pval                 = all_data['sol_Pfvalues'][branch]
                sol_Qval                 = all_data['sol_Qfvalues'][branch]
            elif from_or_to == 't':
                sol_Pval                 = all_data['sol_Ptvalues'][branch]
                sol_Qval                 = all_data['sol_Qtvalues'][branch]

            slack    = coeff_P * sol_Pval + coeff_Q * sol_Qval - z

            if slack > FeasibilityTol:
                log.joint(' WARNING, the limit-envelope cut associated to branch ' + str(branch.count)
                          + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' is violated by the AC solution!\n')
                log.joint(' violation ' + str(slack) + '\n')
                if from_or_to == 'f':
                    log.joint(' values (AC solution) ' + ' Pft ' + str(sol_Pval) + ' Qft ' + str(sol_Qval) + '\n')
                elif from_or_to == 't':
                    log.joint(' values (AC solution) ' + ' Ptf ' + str(sol_Pval) + ' Qtf ' + str(sol_Qval) + '\n')
                #breakexit('check!')
            else:
                log.joint(' AC solution satisfies limit cut at branch ' + str(branch.count) + ' with slack ' + str(slack) + '\n')
            
        
        limit_cuts[(cutid,branch.count)] = (rnd,threshold,from_or_to)
        limit_cuts_info[branch][cutid]   = (rnd,violation,coeff_P,coeff_Q,threshold,cutid,from_or_to)
        
        cutexp = LinExpr()

        if from_or_to == 'f':
            constrname = "limit_cut_" + str(cutid) + "_" + str(branch.count) + "r_" + str(rnd) + "_" + str(f) + "_" + str(t)
            cutexp += coeff_P * Pvar_f[branch] + coeff_Q * Qvar_f[branch]
        elif from_or_to == 't':
            constrname = "limit_cut_" + str(cutid) + "_" + str(branch.count) + "r_" + str(rnd) + "_" + str(t) + "_" + str(f)
            cutexp += coeff_P * Pvar_t[branch] + coeff_Q * Qvar_t[branch]
        else:
            log.joint('  we have a bug\n')
            breakexit('look for bug')
            
        themodel.addConstr(cutexp <= z, name = constrname)

    log.joint('  number limit-envelope cuts added ' + str(most_violated_count) + '\n')
    log.joint('  max error (abs) ' + str(max_error) + ' at branch ' + str(most_violated_branch) + '\n' )

    all_data['most_violated_branch_limit']   = most_violated_branch
    all_data['max_error_limit']              = max_error
    all_data['ID_limit_cuts']               += most_violated_count
    all_data['num_limit_cuts_added']         = most_violated_count
    all_data['num_limit_cuts']              += most_violated_count
    all_data['num_limit_cuts_rnd'][rnd]      = most_violated_count


    if all_data['droplimit']:
        if all_data['addcuts']:
            t0_drop = time.time()
            drop_limit(log,all_data)
            t1_drop = time.time()
            log.joint('  time spent on drop limit ' + str(t1_drop - t0_drop) + '\n')
        elif all_data['round'] >= all_data['cut_age_limit']: ######
            t0_drop = time.time()
            drop_limit(log,all_data)
            t1_drop = time.time()
            log.joint('  time spent on drop limit ' + str(t1_drop - t0_drop) + '\n')


def drop_limit(log,all_data):
    
    #log.joint(' **** drop Limit-envelope cuts ****\n')
    log.joint(' dropping old and slack limit-envelope cuts ...\n')

    themodel                = all_data['themodel']
    branches                = all_data['branches']
    buses                   = all_data['buses']
    IDtoCountmap            = all_data['IDtoCountmap']
    current_rnd             = all_data['round']
    cut_age_limit           = all_data['cut_age_limit']
    limit_cuts              = all_data['limit_cuts']
    limit_cuts_info         = all_data['limit_cuts_info']
    drop_limit              = []
    dropped_limit           = all_data['dropped_limit']
    num_limit_cuts_dropped  = all_data['num_limit_cuts_dropped'] 

    for key in limit_cuts.keys():
        cut     = limit_cuts[key] 
        cut_rnd = cut[0]
        cut_age = current_rnd - cut_rnd

        if cut_age <= cut_age_limit: 
            continue 
        
        cutid         = key[0]
        branchid      = key[1]
        branch        = branches[branchid]
        f             = branch.f
        t             = branch.t
        count_of_f    = IDtoCountmap[f]
        count_of_t    = IDtoCountmap[t]
        cut_threshold = cut[1]
        from_or_to    = cut[2]

        if from_or_to == 'f':
            constrname = "limit_cut_" + str(cutid) + "_" + str(branch.count) + "r_" + str(cut_rnd) + "_" + str(f) + "_" + str(t)
        elif from_or_to == 't':
            constrname = "limit_cut_" + str(cutid) + "_" + str(branch.count) + "r_" + str(cut_rnd) + "_" + str(t) + "_" + str(f)
        else:
            log.joint(' we have a bug\n')
            breakexit('look for bug')

        constr        = themodel.getConstrByName(constrname)
        slack         = - constr.getAttr("slack")

        if ( slack < - cut_threshold ):
            drop_limit.append(key)
            themodel.remove(themodel.getConstrByName(constrname))
            if all_data['loud_cuts']:
                log.joint(' --> removed limit-cut\n')
                log.joint(' the cut ' + str(key)
                          + ' was removed from the model\n')
    
    num_drop_limit = len(drop_limit)
    if num_drop_limit:
        all_data['num_limit_cuts_dropped'] = num_drop_limit
        all_data['num_limit_cuts']        -= num_drop_limit
        all_data['total_limit_dropped']   += num_drop_limit

        dropped_limit.extend(drop_limit)

        for key in drop_limit:
            limit_cuts.pop(key)

        log.joint('  the cuts in drop_limit list were removed from dict limit_cuts\n' )

        ##### drop limit-envelope cuts from limit_cuts_info
        for key in drop_limit:
            cutid    = key[0]
            branchid = key[1]
            cuts_branch = limit_cuts_info[branches[branchid]]
            cuts_branch.pop(cutid)
        log.joint('  cuts in drop_limit list were removed from limit_cuts_info\n') 
    else:
        all_data['num_limit_cuts_dropped'] = 0
        log.joint('  no limit-envelope cuts were dropped this round\n')


def jabr_cuts(log,all_data):
        
    log.joint('\n')
    log.joint(' **** Jabr-cuts ****\n')

    themodel       = all_data['themodel']
    IDtoCountmap   = all_data['IDtoCountmap']
    buses          = all_data['buses']
    branches       = all_data['branches']
    cvar           = all_data['cvar']
    svar           = all_data['svar']
    FeasibilityTol = all_data['tolerance']
    cvalues        = all_data['cvalues']
    svalues        = all_data['svalues']
    
    rnd                 = all_data['round']
    jabr_cuts           = all_data['jabr_cuts']
    jabr_cuts_info      = all_data['jabr_cuts_info']
    num_cuts            = all_data['ID_jabr_cuts']
    num_jabr_cuts_added = 0
    threshold           = all_data['threshold']
    violated            = {}
    violated_count      = 0

    
    all_data['NO_jabrs_violated'] = 0
    
    t0_violation = time.time()

    log.joint(' checking for violations of Jabr inequalities ... \n')

    for branch in branches.values():
        f = branch.f
        t = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        violation = cvalues[branch]*cvalues[branch] + svalues[branch]*svalues[branch] - cvalues[buses[count_of_f]]*cvalues[buses[count_of_t]]
        if violation > threshold:
            violated_count += 1
            violated[branch] = violation

    if violated_count == 0:
        all_data['NO_jabrs_violated'] = 1
        log.joint(' all violations below threshold\n' )
        log.joint(' no more cuts to add for current threshold\n' )
        return None

    t1_violation = time.time()

    log.joint('  number violated Jabrs ' + str(violated_count) + '\n')  
    log.joint('  time spent on violation Jabrs ' + str(t1_violation - t0_violation) + '\n')

    log.joint(' sorting most violated Jabr-envelope cuts ... \n')

    t0_mostviol = time.time()

    num_selected        =  math.ceil(violated_count * all_data['most_violated_fraction_jabr'] )
    most_violated       = dict(sorted(violated.items(), key = lambda x: x[1], reverse = True)[:num_selected])
    most_violated_count = 0

    t1_mostviol = time.time()

    log.joint('  time spent on sorting most violated Jabrs ' + str(t1_mostviol - t0_mostviol) + '\n')

    log.joint(' computing Jabr-envelope cuts ... \n')
    t0_compute  = time.time()

    for branch in most_violated.keys():
        if (most_violated_count == 0):
            most_violated_branch = branch.count
            max_error = most_violated[branch]

        f          = branch.f
        t          = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        cft        = cvalues[branch]
        sft        = svalues[branch]
        cff        = cvalues[buses[count_of_f]]
        ctt        = cvalues[buses[count_of_t]]

        cutnorm   = math.sqrt( (2 * cft)**2 + (2 * sft)**2 + (cff - ctt)**2 )
        coeff_cft = 4 * cft
        coeff_sft = 4 * sft
        coeff_cff = cff - ctt - cutnorm
        coeff_ctt = - (cff - ctt) - cutnorm

        if all_data['parallel_check']:
            if parallel_check(log,all_data,branch,coeff_cft,coeff_sft,coeff_cff,
                              coeff_ctt):
                continue

        #We add the cut
        most_violated_count += 1
        violation            = most_violated[branch]
        cutid                = num_cuts + most_violated_count

        if all_data['loud_cuts']:
            log.joint('  --> new cut\n')
            log.joint('  branch ' + str(branch.count) + ' f ' + str(f) + ' t ' 
                      + str(t) + ' violation ' + str(violation) + ' cut id ' 
                      + str(cutid) + '\n' )
            log.joint('  values ' + ' cft ' + str(cft) + ' sft ' + str(sft) 
                      + ' cff ' + str(cff) + ' ctt ' + str(ctt) + '\n' )
            log.joint('  LHS coeff ' + ' cft ' + str(coeff_cft) + ' sft ' 
                      + str(coeff_sft) + ' cff ' + str(coeff_cff) + ' ctt ' 
                      + str(coeff_ctt) + '\n' )
            log.joint('  cutnorm ' + str(cutnorm) + '\n')

        #sanity check
        if all_data['jabr_validity']:
            
            sol_c     = all_data['sol_cvalues'][branch]
            sol_s     = all_data['sol_svalues'][branch]
            sol_cbusf = all_data['sol_cvalues'][buses[count_of_f]]
            sol_cbust = all_data['sol_cvalues'][buses[count_of_t]]
            slack     = coeff_cft * sol_c + coeff_sft * sol_s + coeff_cff * sol_cbusf + coeff_ctt * sol_cbust

            if slack > FeasibilityTol:
                log.joint('  this cut is not valid!\n')
                log.joint('  violation ' + str(slack) + '\n')
                log.joint('  values (a primal bound)' + ' cft ' + str(sol_c) 
                          + ' sft ' + str(sol_s) + ' cff ' + str(sol_busf) 
                          + ' ctt ' + str(sol_cbust) + '\n' )
                breakexit('check!')
            else:
                log.joint('  valid cut at branch ' + str(branch.count) + ' with slack ' + str(slack) + '\n')
        
        jabr_cuts[(cutid,branch.count)] = (rnd,threshold)
        jabr_cuts_info[branch][cutid]   = (rnd,violation,coeff_cft,coeff_sft,coeff_cff,coeff_ctt,threshold,cutid)
        
        cutexp = LinExpr()
        constrname = "jabr_cut_"+str(cutid)+"_"+str(branch.count)+"r_"+str(rnd)+"_"+str(f)+"_"+str(t)
        cutexp += coeff_cft * cvar[branch] + coeff_sft * svar[branch] + coeff_cff * cvar[buses[count_of_f]] + coeff_ctt * cvar[buses[count_of_t]]

        themodel.addConstr(cutexp <= 0, name = constrname)

    log.joint('  number Jabr-envelope cuts added ' + str(most_violated_count) + '\n')
    log.joint('  max error (abs) ' + str(max_error) + ' at branch ' + str(most_violated_branch) + '\n' )

    t1_compute = time.time()

    log.joint('  time spent on computing Jabrs ' + str(t1_compute - t0_compute) + '\n')

    all_data['most_violated_branch_jabr'] = most_violated_branch
    all_data['max_error_jabr']            = max_error
    all_data['ID_jabr_cuts']             += most_violated_count
    all_data['num_jabr_cuts_added']       = most_violated_count
    all_data['num_jabr_cuts']            += most_violated_count
    all_data['num_jabr_cuts_rnd'][rnd]    = most_violated_count

    if all_data['dropjabr']:
        if all_data['addcuts']:
            t0_drop = time.time()
            drop_jabr(log,all_data)
            t1_drop = time.time()
            log.joint('  time spent on drop Jabrs ' + str(t1_drop - t0_drop) + '\n')
        elif all_data['round'] >= all_data['cut_age_limit']: ######
            t0_drop = time.time()
            drop_jabr(log,all_data)
            t1_drop = time.time()
            log.joint('  time spent on drop Jabrs ' + str(t1_drop - t0_drop) + '\n')

def drop_jabr(log,all_data):
    
    #log.joint('\n')
    #log.joint(' **** drop Jabr-envelope cuts ****\n')
    log.joint(' dropping old and slack Jabr-envelope cuts ...\n')

    themodel               = all_data['themodel']
    branches               = all_data['branches']
    buses                  = all_data['buses']
    IDtoCountmap           = all_data['IDtoCountmap']
    current_rnd            = all_data['round']
    cut_age_limit          = all_data['cut_age_limit']
    jabr_cuts              = all_data['jabr_cuts']
    jabr_cuts_info         = all_data['jabr_cuts_info']
    drop_jabrs             = []
    num_jabr_cuts_dropped  = all_data['num_jabr_cuts_dropped'] 
    dropped_jabrs          = all_data['dropped_jabrs']

    for key in jabr_cuts.keys():
        cut     = jabr_cuts[key] 
        cut_rnd = cut[0]
        cut_age = current_rnd - cut_rnd 

        if cut_age <= cut_age_limit: 
            continue 

        cutid         = key[0]
        branchid      = key[1]
        branch        = branches[branchid]
        f             = branch.f
        t             = branch.t
        count_of_f    = IDtoCountmap[f]
        count_of_t    = IDtoCountmap[t]
        cut_threshold = cut[1]
        
        constrname    = "jabr_cut_"+str(cutid)+"_"+str(branchid)+"r_"+str(cut_rnd)+"_"+str(f)+"_"+str(t)
        constr        = themodel.getConstrByName(constrname)
        slack         = - constr.getAttr("slack")
        
        if ( slack < - cut_threshold ):
            drop_jabrs.append(key)
            themodel.remove(themodel.getConstrByName(constrname))
            if all_data['loud_cuts']:
                log.joint('  --> removed Jabr-cut\n')
                log.joint('  the cut ' + str(key)
                          + ' was removed from the model\n')
    
    num_drop_jabrs = len(drop_jabrs)

    if num_drop_jabrs:
        all_data['num_jabr_cuts_dropped'] = num_drop_jabrs
        all_data['num_jabr_cuts']        -= num_drop_jabrs
        all_data['total_jabr_dropped']   += num_drop_jabrs

        dropped_jabrs.extend(drop_jabrs)

        for key in drop_jabrs:
            jabr_cuts.pop(key)

        log.joint('  the cuts in drop_jabrs list were removed from dict jabr_cuts\n')

        ##### drop Jabr-envelope cuts from jabr_cuts_info
        for key in drop_jabrs:
            cutid    = key[0]
            branchid = key[1]
            cuts_branch = jabr_cuts_info[branches[branchid]]
            cuts_branch.pop(cutid)
        log.joint('  cuts in drop_jabrs list were removed from jabr_cuts_info\n') 
    else:
        all_data['num_jabr_cuts_dropped'] = 0
        log.joint('  no Jabr-envelope cuts were dropped this round\n')


def add_cuts(log,all_data):

    if '_b' in all_data['casename']:
        original_casename = all_data['casename'][:len(all_data['casename']) - 2]
    elif ('_n_5_5' in all_data['casename'] or '_n_0_5' in all_data['casename'] 
          or '_n_1_1' in all_data['casename'] or '_pline' in all_data['casename']):
        original_casename = all_data['casename'][:len(all_data['casename']) - 6]
    else:
        original_casename = all_data['casename']

    filename = '../data/cuts/cuts_' + original_casename + '.txt'
    log.joint(" opening file with cuts " + filename + "\n")

    try:
        thefile = open(filename, "r") 
        lines = thefile.readlines()
    except:
        log.stateandquit(" cannot open file", filename)
        sys.exit("failure")

    numlines  = len(lines)
    theround  = lines[0].split()[3]
    firstline = lines[1].split()
    jabr      = 1
    
    log.joint(' loading cuts from round ' + str(theround) + '\n')

    if firstline[0] == '#Jabr-envelope':
        numjabrcuts = int(firstline[3])
        all_data['addcuts_numjabrcuts'] = numjabrcuts
        log.joint(' number of Jabr-envelope cuts = ' + str(numjabrcuts) + '\n')
    elif firstline[0] == '#i2-envelope':
        jabr      = 0
        i2        = 1
        numi2cuts = int(firstline[3])
        all_data['addcuts_numi2cuts'] = numi2cuts
        log.joint(' no Jabr-envelope cuts to add\n')
        log.joint(' number of i2-envelope cuts = ' + str(numi2cuts) + '\n')
    elif firstline[0] == '#limit-envelope':
        jabr         = 0
        i2           = 0
        numlimitcuts = int(firstline[3])
        all_data['addcuts_numlimitcuts'] = numlimitcuts
        log.joint(' no Jabr nor i2-envelope cuts to add\n')
        log.joint(' number of limit-envelope cuts = ' + str(numlimitcuts) + '\n')
    else:
        log.joint(' no cuts added\n')
        return None

    linenum = 2

    themodel       = all_data['themodel']
    buses          = all_data['buses']
    branches       = all_data['branches']
    IDtoCountmap   = all_data['IDtoCountmap']
    cvar           = all_data['cvar']
    svar           = all_data['svar']
    Pvar_f         = all_data['Pvar_f']
    Qvar_f         = all_data['Qvar_f']
    Pvar_t         = all_data['Pvar_t']
    Qvar_t         = all_data['Qvar_t']    
    FeasibilityTol = all_data['tolerance']

    if all_data['i2']:
        i2var_f = all_data['i2var_f']

    while linenum < numlines: 
        thisline = lines[linenum].split()
        if thisline[0] == '#i2-envelope' and jabr:
            numi2cuts = int(thisline[3])
            all_data['addcuts_numi2cuts'] = numi2cuts
            log.joint(' number of i2-envelope cuts = ' + str(numi2cuts) + '\n')
            linenum += 1
            jabr = 0
            i2   = 1
            continue

        elif thisline[0] == '#limit-envelope' and i2:
            numlimitcuts = int(thisline[3])
            all_data['addcuts_numlimitcuts'] = numlimitcuts
            log.joint(' number of limit-envelope cuts = ' + str(numlimitcuts) 
                      + '\n')
            linenum += 1
            i2       = 0
            continue

        elif jabr:
            branchid  = int(thisline[1])
            f         = int(thisline[3])
            t         = int(thisline[5])
            cutid     = int(thisline[7])
            rnd       = int(thisline[9])
            coeff_cft = float(thisline[15])
            coeff_sft = float(thisline[17])
            coeff_cff = float(thisline[19])
            coeff_ctt = float(thisline[21])


            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' 
                          + str(branchid) + ' f ' + str(f) + ' t ' + str(t) 
                          + ' was turned OFF\n')
                linenum += 1
                continue

            branch     = branches[branchid]
            count_of_f = IDtoCountmap[f] 
            count_of_t = IDtoCountmap[t]

            if ( (branchid != branch.count) or f != (branch.f) 
                 or (t != branch.t) ):

                log.joint(' branchid ' + str(branchid) + ' branch.count ' 
                          + str(branch.count) + ' f ' + str(f) + ' branch.f ' 
                          + str(branch.f) + ' t ' + str(t) + ' branch.t ' 
                          + str(branch.t) + '\n')
                breakexit('bug')

            if all_data['loud_cuts']:
                log.joint(' --> new Jabr-envelope cut\n')
                log.joint(' branch ' + str(branchid) + ' f ' + str(f) + ' t ' 
                          + str(t) + ' cutid ' + str(cutid) + '\n' )
                log.joint(' LHS coeff ' + ' cft ' + str(coeff_cft) + ' sft ' 
                          + str(coeff_sft) + ' cff ' + str(coeff_cff) 
                          + ' ctt ' + str(coeff_ctt) + '\n' )
                            
            if all_data['jabr_validity']:
                sol_c           = all_data['sol_cvalues'][branch]
                sol_s           = all_data['sol_svalues'][branch]
                sol_cbusf       = all_data['sol_cvalues'][buses[count_of_f]]
                sol_cbust       = all_data['sol_cvalues'][buses[count_of_t]]
                sol_violation   = coeff_cft * sol_c + coeff_sft * sol_s + coeff_cff * sol_cbusf + coeff_ctt * sol_cbust
                sol_relviolation   = sol_violation / ( ( coeff_cft**2 + coeff_sft**2 + coeff_cff**2 + coeff_ctt**2 )**0.5 ) 

                if sol_relviolation > FeasibilityTol:
                    log.joint(' WARNING, the Jabr-envelope cut associated to branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' is violated by the AC solution!\n')
                    log.joint(' violation ' + str(sol_violation) + '\n')
                    log.joint(' relative violation ' + str(sol_relviolation) + '\n')
                    log.joint(' values (AC solution)' + ' cft ' + str(sol_c) + ' sft ' + str(sol_s) + ' cff ' + str(sol_busf) + ' ctt ' + str(sol_cbust) + '\n' )
                    breakexit('check!')
                else:
                    log.joint(' AC solution satisfies Jabr inequality at branch ' + str(branch.count) + ' with slack ' + str(sol_relviolation) + '\n')


            cutexp = LinExpr()
            constrname = "jabr_cut_"+str(cutid)+"_"+str(branchid)+"r_"+str(rnd)+"_"+str(f)+"_"+str(t)
            cutexp += coeff_cft * cvar[branch] + coeff_sft * svar[branch] + coeff_cff * cvar[buses[count_of_f]] + coeff_ctt * cvar[buses[count_of_t]]

            themodel.addConstr(cutexp <= 0, name = constrname)            
            linenum += 1

        elif (jabr == 0) and i2:

            branchid   = int(thisline[1])
            f          = int(thisline[3])
            t          = int(thisline[5])
            cutid      = int(thisline[7])
            coeff_Pft  = float(thisline[15])
            coeff_Qft  = float(thisline[17])
            coeff_cff  = float(thisline[19])
            coeff_i2ft = float(thisline[21])

            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' 
                          + str(branchid) + ' f ' + str(f) + ' t ' + str(t) 
                          + 'was turned OFF\n')
                linenum += 1
                continue

            branch     = branches[branchid]
            count_of_f = IDtoCountmap[f] 
            count_of_t = IDtoCountmap[t]

            if ( (branchid != branch.count) or f != (branch.f) 
                 or (t != branch.t) ):
                breakexit('there might be bug')

            if all_data['loud_cuts']:
                log.joint(' --> new i2-envelope cut\n')
                log.joint(' branch ' + str(branchid) + ' f ' + str(f) + ' t ' 
                          + str(t) + ' cutid ' + str(cutid) + '\n' )
                log.joint(' LHS coeff ' + ' Pft ' + str(coeff_Pft) + ' Qft ' 
                          + str(coeff_Qft) + ' cff ' + str(coeff_cff) 
                          + ' i2ft ' + str(coeff_i2ft) + '\n' )
            
            if all_data['i2_validity']:
                sol_Pf           = all_data['sol_Pfvalues'][branch]
                sol_Qf           = all_data['sol_Qfvalues'][branch]
                sol_c            = all_data['sol_cvalues'][branch]
                sol_s            = all_data['sol_svalues'][branch]
                sol_cbusf        = all_data['sol_cvalues'][buses[count_of_f]]
                sol_cbust        = all_data['sol_cvalues'][buses[count_of_t]]
                sol_i2f          = computei2value(log,all_data,branch,sol_c,sol_s,sol_cbusf,sol_cbust)
                sol_violation    = coeff_Pft * sol_Pf + coeff_Qft * sol_Qf + coeff_cff * sol_cbusf + coeff_i2ft * sol_i2f
                sol_relviolation = sol_violation / ( ( coeff_Pft**2 + coeff_Qft**2 + coeff_cff**2 + coeff_i2ft**2 )**0.5 ) 

                if sol_relviolation > FeasibilityTol:
                    log.joint(' WARNING, the i2-envelope cut associated to branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' is violated by the AC solution!\n')
                    log.joint(' violation ' + str(sol_violation) + '\n')
                    log.joint(' relative violation ' + str(sol_relviolation) + '\n')
                    log.joint(' values (AC solution) ' + ' Pft ' + str(sol_Pf) + ' Qft ' + str(sol_Qf) + ' cff ' + str(sol_cbusf) + ' i2ft ' + str(sol_i2f) + '\n' )
                    breakexit('check!')
                else:
                    log.joint(' AC solution satisfies i2 cut at branch ' + str(branch.count) + ' with slack ' + str(sol_relviolation) + '\n')


            cutexp     = LinExpr()
            constrname = "i2_cut_"+str(cutid)+"_"+str(branch.count)+"r_"+str(rnd)+"_"+str(f)+"_"+str(t)
            cutexp    += coeff_Pft * Pvar_f[branch] + coeff_Qft * Qvar_f[branch] + coeff_cff * cvar[buses[count_of_f]] + coeff_i2ft * i2var_f[branch]

            themodel.addConstr(cutexp <= 0, name = constrname)
            linenum += 1


        elif (jabr == 0) and (i2 == 0):
            branchid   = int(thisline[1])
            f          = int(thisline[3])
            t          = int(thisline[5])
            cutid      = int(thisline[7])
            if thisline[14] == 'Pft':
                from_or_to = 'f'
            elif thisline[14] == 'Ptf':
                from_or_to = 't'
            else:
                log.joint(' look for a bug\n')
                breakexit('bug')

            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' 
                          + str(branchid) + ' f ' + str(f) + ' t ' + str(t) 
                          + 'was turned OFF\n')
                linenum += 1
                continue

            coeff_P    = float(thisline[15])
            coeff_Q    = float(thisline[17])

            branch     = branches[branchid]
            count_of_f = IDtoCountmap[f] 
            count_of_t = IDtoCountmap[t]

            if ( (branchid != branch.count) or f != (branch.f) 
                 or (t != branch.t) ):
                breakexit('there might be bug')

            if all_data['loud_cuts']:
                log.joint(' --> new limit-envelope cut\n')
                log.joint(' branch ' + str(branchid) + ' f ' + str(f) + ' t ' 
                          + str(t) + ' cutid ' + str(cutid) + '\n' )
                if from_or_to == 'f':
                    log.joint(' LHS coeff ' + ' Pft ' + str(coeff_P) 
                              + ' Qft ' + str(coeff_Q) + '\n')
                elif from_or_to == 't':
                    log.joint(' LHS coeff ' + ' Ptf ' + str(coeff_P) 
                              + ' Qtf ' + str(coeff_Q) + '\n')
            
            
            if all_data['limit_validity']:
                if from_or_to == 'f':
                    sol_P           = all_data['sol_Pfvalues'][branch]
                    sol_Q           = all_data['sol_Qfvalues'][branch]
                elif from_or_to == 't':
                    sol_P           = all_data['sol_Ptvalues'][branch]
                    sol_Q           = all_data['sol_Qtvalues'][branch]
                sol_violation    = coeff_P * sol_P + coeff_Q * sol_Q - 1

                if sol_violation > FeasibilityTol:
                    log.joint(' WARNING, the limit-envelope cut associated to branch ' + str(branch.count)
                              + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' is violated by the AC solution!\n')
                    log.joint(' violation ' + str(sol_violation) + '\n')
                    if from_or_to == 'f':
                        log.joint(' values (AC solution) ' + ' Pft ' + str(sol_P) + ' Qft ' + str(sol_Q) + '\n')
                    elif from_or_to == 't':
                        log.joint(' values (AC solution) ' + ' Ptf ' + str(sol_P) + ' Qtf ' + str(sol_Q) + '\n')
                    #breakexit('check!')
                else:
                    log.joint(' AC solution satisfies limit cut at branch ' + str(branch.count) + ' with slack ' + str(sol_violation) + '\n')


            cutexp     = LinExpr()
            if from_or_to == 'f':
                constrname = "limit_cut_"+str(cutid)+"_"+str(branch.count)+"r_"+str(rnd)+"_"+str(f)+"_"+str(t)
                cutexp    += coeff_P * Pvar_f[branch] + coeff_Q * Qvar_f[branch]
            elif from_or_to == 't':
                constrname = "limit_cut_"+str(cutid)+"_"+str(branch.count)+"r_"+str(rnd)+"_"+str(f)+"_"+str(t)
                cutexp    += coeff_P * Pvar_t[branch] + coeff_Q * Qvar_t[branch]

            themodel.addConstr(cutexp <= 1, name = constrname)
            linenum += 1

    all_data['num_jabr_cuts']  = numjabrcuts
    all_data['num_i2_cuts']    = numi2cuts
    all_data['num_limit_cuts'] = numlimitcuts
    all_data['ID_jabr_cuts']   = numjabrcuts
    all_data['ID_i2_cuts']     = numi2cuts
    all_data['ID_limit_cuts']  = numlimitcuts


def write_cuts(log,all_data):

    filename = 'cuts_' + all_data['casename'] + '.txt' 
    log.joint(" opening file with cuts " + filename + "\n")

    try:
        thefile = open(filename, "w+")
        lines = thefile.readlines()
    except:
        log.stateandquit(" cannot open file", filename)
        sys.exit("failure")
    

    jabr_cuts = all_data['jabr_cuts_info']
    rnd       = all_data['round'] 

    log.joint(' writing down to ' + filename + ' jabr-envelope cuts in round ' + str(rnd) + '\n')

    thefile.write('current round = ' + str(all_data['round']) + '\n')
    thefile.write('#Jabr-envelope cuts = ' + str(all_data['num_jabr_cuts']) + '\n')
    for branch in jabr_cuts.keys():
        branchcount = branch.count
        f           = branch.f
        t           = branch.t
        branch_cuts = jabr_cuts[branch]
        
        for cutid in branch_cuts.keys():
            cut = branch_cuts[cutid]
            if len(cut) == 0:
                continue
            cut_info = 'branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' cutid ' + str(cutid) + ' round ' + str(cut[0]) + ' violation ' + str(cut[1]) + ' threshold ' + str(cut[6]) + ' cft ' + str(cut[2]) + ' sft ' + str(cut[3]) + ' cff ' + str(cut[4]) + ' ctt ' + str(cut[5]) + '\n'
            thefile.write(cut_info)
        
    log.joint(' done writing down jabr-envelope cuts\n')

    if all_data['i2cuts']:
        log.joint(' writing down to ' + filename + ' i2-envelope cuts in round ' + str(rnd) + '\n')
        thefile.write('#i2-envelope cuts = ' + str(all_data['num_i2_cuts']) + '\n')
        
        i2_cuts = all_data['i2_cuts_info']
        for branch in i2_cuts.keys(): #BUG found
            branchcount = branch.count
            f           = branch.f
            t           = branch.t
            branch_cuts = i2_cuts[branch]
        
            for cutid in branch_cuts.keys():
                cut = branch_cuts[cutid]
                if len(cut) == 0:
                    continue
                cut_info = 'branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' cutid ' + str(cutid) + ' round ' + str(cut[0]) + ' violation ' + str(cut[1]) + ' threshold ' + str(cut[6]) + ' Pft ' + str(cut[2]) + ' Qft ' + str(cut[3]) + ' cff ' + str(cut[4]) + ' i2ft ' + str(cut[5]) + '\n'
                thefile.write(cut_info)        

        log.joint(' done writing down i2-envelope cuts\n')


    if all_data['limitcuts']:
        log.joint(' writing down to ' + filename + ' limit-envelope cuts in round ' + str(rnd) + '\n')
        thefile.write('#limit-envelope cuts = ' + str(all_data['num_limit_cuts']) + '\n')
        
        limit_cuts = all_data['limit_cuts_info']
        for branch in limit_cuts.keys():
            branchcount = branch.count
            f           = branch.f    #(rnd,violation,coeff_P,coeff_Q,threshold,cutid,from_or_to)
            t           = branch.t
            branch_cuts = limit_cuts[branch]
        
            for cutid in branch_cuts.keys():
                cut        = branch_cuts[cutid]

                if len(cut) == 0:
                    continue

                from_or_to = cut[6]
                if from_or_to == 'f':
                    cut_info = 'branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' cutid ' + str(cutid) + ' round ' + str(cut[0]) + ' violation ' + str(cut[1]) + ' threshold ' + str(cut[4]) + ' Pft ' + str(cut[2]) + ' Qft ' + str(cut[3]) + '\n'
                elif from_or_to == 't':
                    cut_info = 'branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' cutid ' + str(cutid) + ' round ' + str(cut[0]) + ' violation ' + str(cut[1]) + ' threshold ' + str(cut[4]) + ' Ptf ' + str(cut[2]) + ' Qtf ' + str(cut[3]) + '\n'
                else:
                    log.joint(' we have a bug\n')
                    breakexit('look for the bug')

                thefile.write(cut_info)        

        log.joint(' done writing down limit-envelope cuts\n')


    thefile.close()

def parallel_check(log,all_data,branch,coeff_cft,coeff_sft,coeff_cff,coeff_ctt): 

    threshold         = all_data['threshold']
    threshold_dotprod = all_data['threshold_dotprod']    
    jabr_cuts_info    = all_data['jabr_cuts_info']
    cuts_branch       = jabr_cuts_info[branch] 
  
    if all_data['loud_cuts']:
        log.joint('\n -- parallel check wrt previous Jabr-envelope cuts at branch ' + str(branch.count) +'\n')

    if len(cuts_branch) == 0:
        if all_data['loud_cuts']:
            log.joint(' first Jabr-envelope cut, we add it\n')
        return 0 
    else:
        v = compute_normal(coeff_cft,coeff_sft,coeff_cff,coeff_ctt)
        if all_data['loud_cuts']:
            log.joint(' LHS coeffs potential cut ' + ' cft ' + str(coeff_cft) 
                      + ' sft ' + str(coeff_sft) + ' cff ' + str(coeff_cff) 
                      + ' ctt ' + str(coeff_ctt) + '\n' )

        for cut in cuts_branch.values():
            cutid         = cut[7]
            cut_coeff_cft = cut[2]
            cut_coeff_sft = cut[3]
            cut_coeff_cff = cut[4]
            cut_coeff_ctt = cut[5]
            w = compute_normal(cut_coeff_cft,cut_coeff_sft,cut_coeff_cff,
                               cut_coeff_ctt)
            dotprod = np.dot(v,w)
            angle = np.arccos(dotprod)
            angle_deg = angle * 180 / np.pi
            if all_data['loud_cuts']:
                log.joint(' LHS coeffs of cutid ' + str(cutid) + ' cft ' 
                          + str(cut_coeff_cft) + ' sft ' + str(cut_coeff_sft) 
                          + ' cff ' + str(cut_coeff_cff) + ' ctt ' 
                          + str(cut_coeff_ctt) + '\n' )
                log.joint(' angle (rad) ' + str(angle) + ' angle (deg) ' 
                          + str(angle_deg) +  ' dot-product ' + str(dotprod) 
                          + '\n')
            
            if dotprod > 1 - threshold_dotprod:
                if all_data['loud_cuts']:
                    log.joint(' parallel cut, should not be added\n')
                return 1
            else:
                if all_data['loud_cuts']:
                    log.joint(' cut should be added\n')
                return 0
        
     
def parallel_check_i2(log,all_data,branch,coeff_Pft,coeff_Qft,coeff_cff,coeff_i2ft):

    threshold         = all_data['threshold']
    threshold_dotprod = all_data['threshold_dotprod']    
    i2_cuts_info      = all_data['i2_cuts_info']      
    cuts_branch       = i2_cuts_info[branch]          

    if all_data['loud_cuts']:
        log.joint('\n -- parallel check wrt previous i2-envelope cuts at branch ' + str(branch.count) +'\n')

    if len(cuts_branch) == 0:
        if all_data['loud_cuts']:
            log.joint(' first i2-envelope cut, we add it\n')
        return 0 
    else:
        v = compute_normal(coeff_Pft,coeff_Qft,coeff_cff,coeff_i2ft)
        if all_data['loud_cuts']:
            log.joint(' LHS coeffs potential cut ' + ' Pft ' + str(coeff_Pft) 
                      + ' Qft ' + str(coeff_Qft) + ' cff ' + str(coeff_cff) 
                      + ' i2ft ' + str(coeff_i2ft) + '\n' )

        for cut in cuts_branch.values():
            cutid      = cut[7]
            cut_coeff_Pft  = cut[2]
            cut_coeff_Qft  = cut[3]
            cut_coeff_cff  = cut[4]
            cut_coeff_i2ft = cut[5]
            w = compute_normal(cut_coeff_Pft,cut_coeff_Qft,cut_coeff_cff,
                               coeff_i2ft)
            dotprod = np.dot(v,w)
            angle = np.arccos(dotprod)
            angle_deg = angle * 180 / np.pi
            if all_data['loud_cuts']:
                log.joint(' LHS coeffs of cutid ' + str(cutid) + ' cft ' 
                          + str(cut_coeff_Pft) + ' sft ' + str(cut_coeff_Qft) 
                          + ' cff ' + str(cut_coeff_cff) + ' ctt ' 
                          + str(cut_coeff_i2ft) + '\n' )
                log.joint(' angle (rad) ' + str(angle) + ' angle (deg) ' 
                          + str(angle_deg) +  ' dot-product ' + str(dotprod) 
                          + '\n')
            if dotprod > 1 - threshold_dotprod:
                if all_data['loud_cuts']:
                    log.joint(' parallel cut, should not be added\n')
                return 1
            else:
                if all_data['loud_cuts']:
                    log.joint(' cut should be added\n')
                return 0


def parallel_check_limit(log,all_data,branch,coeff_P,coeff_Q,from_or_to):

    threshold         = all_data['threshold']
    threshold_dotprod = all_data['threshold_dotprod']    
    limit_cuts_info   = all_data['limit_cuts_info']
    cuts_branch       = limit_cuts_info[branch]          

    if all_data['loud_cuts']:
        log.joint('\n -- parallel check wrt previous limit-envelope cuts at branch ' + str(branch.count) +'\n')

    if len(cuts_branch) == 0:
        if all_data['loud_cuts']:
            log.joint(' first limit-envelope cut, we add it\n')
        return 0 
    else:        
        v = compute_normal(coeff_P,coeff_Q)
        if all_data['loud_cuts']:
            if from_or_to == 'f':
                log.joint(' LHS coeffs potential cut ' + ' Pft ' 
                          + str(coeff_P) + ' Qft ' + str(coeff_Q) + '\n')
            elif from_or_to == 't':
                log.joint(' LHS coeffs potential cut ' + ' Ptf ' 
                          + str(coeff_P) + ' Qtf ' + str(coeff_Q) + '\n' )

        for cut in cuts_branch.values(): 
            cutid          = cut[5]
            cut_coeff_P    = cut[2]
            cut_coeff_Q    = cut[3]
            cut_from_or_to = cut[6]

            if from_or_to != cut_from_or_to: #we check whether potential cut has the same sense (from/to) as the incumbent cut
                continue

            w         = compute_normal(cut_coeff_P,cut_coeff_Q)
            dotprod   = np.dot(v,w)
            angle     = np.arccos(dotprod)
            angle_deg = angle * 180 / np.pi
            
            if all_data['loud_cuts']:
                if cut_from_or_to == 'f':
                    log.joint(' LHS coeffs of cutid ' + str(cutid) + ' Pft ' 
                              + str(cut_coeff_P) + ' Qft ' + str(cut_coeff_Q) 
                              + '\n')
                elif cut_from_or_to == 't':
                    log.joint(' LHS coeffs of cutid ' + str(cutid) + ' Ptf ' 
                              + str(cut_coeff_P) + ' Qtf ' + str(cut_coeff_Q) 
                              + '\n')
                log.joint(' angle (rad) ' + str(angle) + ' angle (deg) ' 
                          + str(angle_deg) +  ' dot-product ' + str(dotprod) 
                          + '\n')

            if dotprod > 1 - threshold_dotprod:
                if all_data['loud_cuts']:
                    log.joint(' parallel cut, should not be added\n')
                return 1
            else:
                if all_data['loud_cuts']:
                    log.joint(' cut should be added\n')
                return 0
            
        




def computei2value(log,all_data,branch,sol_c,sol_s,sol_cbusf,sol_cbust):
    
    ratio  = branch.ratio
    y      = branch.y
    g      = y.real
    b      = y.imag
    bshunt = branch.bc
    angle  = branch.angle_rad
                                                                                
    #from-to
    i2f  = 0
    i2f += (g*g + b*b)/(ratio*ratio) * ( (sol_cbusf/(ratio*ratio)) + sol_cbust - (2/ratio) * ( sol_c * math.cos(angle) + sol_s * math.sin(angle) ) )
    i2f += b*bshunt/(ratio**3) * ( (sol_cbusf/ratio) - (sol_c * math.cos(angle) + sol_s * math.sin(angle) )) 
    i2f += g*bshunt/(ratio**3) * ( sol_s * math.cos(angle) - sol_c * math.sin(angle) )
    i2f += (bshunt*bshunt*sol_cbusf/(4*(ratio**4)) )
        
    return i2f



def add_cuts_ws(log,all_data):

    if '_b' in all_data['casename']:
        original_casename = all_data['casename'][:len(all_data['casename']) - 2]
    elif ('_n_5_5' in all_data['casename'] or '_n_0_5' in all_data['casename'] 
          or '_n_1_1' in all_data['casename'] or '_pline' in all_data['casename']):
        original_casename = all_data['casename'][:len(all_data['casename']) - 6]        
    #elif '_line5' in all_data['casename']:
    #    original_casename = all_data['casename'][:len(all_data['casename']) - 6]
    else:
        original_casename = all_data['casename']

    filename = '../data/cuts/cuts_' + original_casename + '.txt'
    log.joint(" opening file with cuts " + filename + "\n")


    try:
        thefile = open(filename, "r") 
        lines = thefile.readlines()
    except:
        log.stateandquit(" cannot open file", filename)
        sys.exit("failure")

    numlines  = len(lines)
    theround  = lines[0].split()[3]
    firstline = lines[1].split()
    jabr      = 1
    
    log.joint(' loading cuts from round ' + str(theround) + '\n')

    if firstline[0] == '#Jabr-envelope':
        numjabrcuts = int(firstline[3])
        all_data['addcuts_numjabrcuts'] = numjabrcuts
        log.joint(' number of Jabr-envelope cuts = ' + str(numjabrcuts) + '\n')
    elif firstline[0] == '#i2-envelope':
        jabr      = 0
        i2        = 1
        numi2cuts = int(firstline[3])
        all_data['addcuts_numi2cuts'] = numi2cuts
        log.joint(' no Jabr-envelope cuts to add\n')
        log.joint(' number of i2-envelope cuts = ' + str(numi2cuts) + '\n')
    elif firstline[0] == '#limit-envelope':
        jabr         = 0
        i2           = 0
        numlimitcuts = int(firstline[3])
        all_data['addcuts_numlimitcuts'] = numlimitcuts
        log.joint(' no Jabr nor i2-envelope cuts to add\n')
        log.joint(' number of limit-envelope cuts = ' + str(numlimitcuts) + '\n')
    else:
        log.joint(' no cuts added\n')
        return None

    linenum = 2

    themodel       = all_data['themodel']
    buses          = all_data['buses']
    branches       = all_data['branches']
    IDtoCountmap   = all_data['IDtoCountmap']
    cvar           = all_data['cvar']
    svar           = all_data['svar']
    Pvar_f         = all_data['Pvar_f']
    Qvar_f         = all_data['Qvar_f']
    Pvar_t         = all_data['Pvar_t']
    Qvar_t         = all_data['Qvar_t']    
    FeasibilityTol = all_data['tolerance']

    if all_data['i2']:
        i2var_f = all_data['i2var_f']

    #load dictionaries
    
    jabr_cuts       = {} #(rnd,threshold)
    jabr_cuts_info  = {} #jabr_cuts_info[branch][cutid]   = (rnd,violation,coeff_cft,coeff_sft,coeff_cff,coeff_ctt,threshold,cutid)
    i2_cuts         = {} #(rnd,threshold)
    i2_cuts_info    = {} #(rnd,violation,coeff_Pft,coeff_Qft,coeff_cff,coeff_i2ft,threshold,cutid)
    limit_cuts      = {} #(rnd,threshold,from_or_to)
    limit_cuts_info = {} #(rnd,violation,coeff_P,coeff_Q,threshold,cutid,from_or_to)
    
    for branch in branches.values():
        jabr_cuts_info[branch]  = {}
        i2_cuts_info[branch]    = {}
        limit_cuts_info[branch] = {}
    
    cutid_jabr  = 0
    cutid_i2    = 0
    cutid_limit = 0
        
    while linenum < numlines: 
        thisline = lines[linenum].split()
        if thisline[0] == '#i2-envelope' and jabr:
            numi2cuts = int(thisline[3])
            all_data['addcuts_numi2cuts'] = numi2cuts
            log.joint(' number of i2-envelope cuts = ' + str(numi2cuts) + '\n')
            linenum += 1
            jabr = 0
            i2   = 1
            continue

        elif thisline[0] == '#limit-envelope' and i2:
            numlimitcuts = int(thisline[3])
            all_data['addcuts_numlimitcuts'] = numlimitcuts
            log.joint(' number of limit-envelope cuts = ' + str(numlimitcuts) + '\n')
            linenum += 1
            i2       = 0
            continue

        elif jabr:
            branchid  = int(thisline[1])
            f         = int(thisline[3])
            t         = int(thisline[5])
            #cutid     = int(thisline[7])
            rnd       = int(thisline[9])
            violation = float(thisline[11])
            threshold = float(thisline[13])
            coeff_cft = float(thisline[15])
            coeff_sft = float(thisline[17])
            coeff_cff = float(thisline[19])
            coeff_ctt = float(thisline[21])


            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' + str(branchid) 
                          + ' f ' + str(f) + ' t ' + str(t) 
                          + ' was turned OFF\n')
                linenum += 1
                continue

            branch     = branches[branchid]
            count_of_f = IDtoCountmap[f] 
            count_of_t = IDtoCountmap[t]

            if (branchid != branch.count) or f != (branch.f) or (t != branch.t):
                log.joint(' branchid ' + str(branchid) + ' branch.count ' 
                          + str(branch.count) + ' f ' + str(f) + ' branch.f ' 
                          + str(branch.f) + ' t ' + str(t) + ' branch.t ' 
                          + str(branch.t) + '\n')
                breakexit('bug')

            cutid_jabr += 1

            jabr_cuts[(cutid_jabr,branchid)]   = (rnd,threshold)
            jabr_cuts_info[branch][cutid_jabr] = (rnd,violation,coeff_cft,coeff_sft,coeff_cff,coeff_ctt,threshold,cutid_jabr)

            if all_data['loud_cuts']:
                log.joint(' --> new Jabr-envelope cut\n')
                log.joint(' branch ' + str(branchid) + ' f ' + str(f) + ' t ' 
                          + str(t) + ' cutid ' + str(cutid_jabr) + '\n' )
                log.joint(' LHS coeff ' + ' cft ' + str(coeff_cft) + ' sft ' 
                          + str(coeff_sft) + ' cff ' + str(coeff_cff) 
                          + ' ctt ' + str(coeff_ctt) + '\n' )
                            
            if all_data['jabr_validity']:
                sol_c           = all_data['sol_cvalues'][branch]
                sol_s           = all_data['sol_svalues'][branch]
                sol_cbusf       = all_data['sol_cvalues'][buses[count_of_f]]
                sol_cbust       = all_data['sol_cvalues'][buses[count_of_t]]
                sol_violation   = coeff_cft * sol_c + coeff_sft * sol_s + coeff_cff * sol_cbusf + coeff_ctt * sol_cbust
                sol_relviolation   = sol_violation / ( ( coeff_cft**2 + coeff_sft**2 + coeff_cff**2 + coeff_ctt**2 )**0.5 ) 

                if sol_relviolation > FeasibilityTol:
                    log.joint(' WARNING, the Jabr-envelope cut associated to branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' is violated by the AC solution!\n')
                    log.joint(' violation ' + str(sol_violation) + '\n')
                    log.joint(' relative violation ' + str(sol_relviolation) + '\n')
                    log.joint(' values (AC solution)' + ' cft ' + str(sol_c) + ' sft ' + str(sol_s) + ' cff ' + str(sol_busf) + ' ctt ' + str(sol_cbust) + '\n' )
                    #breakexit('check!')
                else:
                    log.joint(' AC solution satisfies Jabr inequality at branch ' + str(branch.count) + ' with slack ' + str(sol_relviolation) + '\n')


            cutexp = LinExpr()
            constrname = "jabr_cut_"+str(cutid_jabr)+"_"+str(branchid)+"r_"+str(rnd)+"_"+str(f)+"_"+str(t)
            cutexp += coeff_cft * cvar[branch] + coeff_sft * svar[branch] + coeff_cff * cvar[buses[count_of_f]] + coeff_ctt * cvar[buses[count_of_t]]

            themodel.addConstr(cutexp <= 0, name = constrname)            
            linenum += 1

        elif (jabr == 0) and i2:
            #log.joint(' here we have the bug, printing thisline[1] ' + thisline[1] + '\n')
            branchid   = int(thisline[1])
            f          = int(thisline[3])
            t          = int(thisline[5])
            #cutid      = int(thisline[7])
            rnd        = int(thisline[9])
            violation  = float(thisline[11])
            threshold  = float(thisline[13])
            coeff_Pft  = float(thisline[15])
            coeff_Qft  = float(thisline[17])
            coeff_cff  = float(thisline[19])
            coeff_i2ft = float(thisline[21])

            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + 'was turned OFF\n')
                linenum += 1
                continue

            branch     = branches[branchid]
            count_of_f = IDtoCountmap[f] 
            count_of_t = IDtoCountmap[t]

            if (branchid != branch.count) or f != (branch.f) or (t != branch.t):
                breakexit('there might be bug')


            cutid_i2 += 1 # updating i2-cut counter

            # updating dictionaries
            
            i2_cuts[(cutid_i2,branchid)]   = (rnd,threshold)
            i2_cuts_info[branch][cutid_i2] = (rnd,violation,coeff_Pft,coeff_Qft,coeff_cff,coeff_i2ft,threshold,cutid_i2)
            

            if all_data['loud_cuts']:
                log.joint(' --> new i2-envelope cut\n')
                log.joint(' branch ' + str(branchid) + ' f ' + str(f) + ' t ' 
                          + str(t) + ' cutid ' + str(cutid_i2) + '\n' )
                log.joint(' LHS coeff ' + ' Pft ' + str(coeff_Pft) + ' Qft ' 
                          + str(coeff_Qft) + ' cff ' + str(coeff_cff) 
                          + ' i2ft ' + str(coeff_i2ft) + '\n' )
            
            if all_data['i2_validity']:
                sol_Pf           = all_data['sol_Pfvalues'][branch]
                sol_Qf           = all_data['sol_Qfvalues'][branch]
                sol_c            = all_data['sol_cvalues'][branch]
                sol_s            = all_data['sol_svalues'][branch]
                sol_cbusf        = all_data['sol_cvalues'][buses[count_of_f]]
                sol_cbust        = all_data['sol_cvalues'][buses[count_of_t]]
                sol_i2f          = computei2value(log,all_data,branch,sol_c,sol_s,sol_cbusf,sol_cbust)
                sol_violation    = coeff_Pft * sol_Pf + coeff_Qft * sol_Qf + coeff_cff * sol_cbusf + coeff_i2ft * sol_i2f
                sol_relviolation = sol_violation / ( ( coeff_Pft**2 + coeff_Qft**2 + coeff_cff**2 + coeff_i2ft**2 )**0.5 ) 

                if sol_relviolation > FeasibilityTol:
                    log.joint(' WARNING, the i2-envelope cut associated to branch ' 
                              + str(branchid) + ' f ' + str(f) + ' t ' + str(t) 
                              + ' is violated by the AC solution!\n')
                    log.joint(' violation ' + str(sol_violation) + '\n')
                    log.joint(' relative violation ' + str(sol_relviolation) + '\n')
                    log.joint(' values (AC solution) ' + ' Pft ' + str(sol_Pf) 
                              + ' Qft ' + str(sol_Qf) + ' cff ' + str(sol_cbusf) 
                              + ' i2ft ' + str(sol_i2f) + '\n' )
                    #breakexit('check!')
                else:
                    log.joint(' AC solution satisfies i2 inequality at branch ' + str(branchid) + ' with slack ' + str(sol_relviolation) + '\n')

            cutexp     = LinExpr()
            constrname = "i2_cut_"+str(cutid_i2)+"_"+str(branchid)+"r_"+str(rnd)+"_"+str(f)+"_"+str(t)
            cutexp    += coeff_Pft * Pvar_f[branch] + coeff_Qft * Qvar_f[branch] + coeff_cff * cvar[buses[count_of_f]] + coeff_i2ft * i2var_f[branch]

            themodel.addConstr(cutexp <= 0, name = constrname)
            linenum += 1


        elif (jabr == 0) and (i2 == 0):
            branchid   = int(thisline[1])
            f          = int(thisline[3])
            t          = int(thisline[5])
            cutid      = int(thisline[7])
            if thisline[14] == 'Pft':
                from_or_to = 'f'
            elif thisline[14] == 'Ptf':
                from_or_to = 't'
            else:
                log.joint(' look for a bug\n')
                breakexit('bug')

            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + 'was turned OFF\n')
                linenum += 1
                continue

            rnd        = int(thisline[9])
            violation  = float(thisline[11])
            threshold  = float(thisline[13])
            coeff_P    = float(thisline[15])
            coeff_Q    = float(thisline[17])

            branch     = branches[branchid]
            count_of_f = IDtoCountmap[f] 
            count_of_t = IDtoCountmap[t]

            if (branchid != branch.count) or f != (branch.f) or (t != branch.t):
                breakexit('there might be bug')

            cutid_limit += 1 # updating limit-cut counter
            
            # Updating dictionaries
            
            limit_cuts[(cutid_limit,branchid)]   = (rnd,threshold,from_or_to)
            limit_cuts_info[branch][cutid_limit] = (rnd,violation,coeff_P,coeff_Q,threshold,cutid_limit,from_or_to)

            if all_data['loud_cuts']:
                log.joint(' --> new limit-envelope cut\n')
                log.joint(' branch ' + str(branchid) + ' f ' + str(f) + ' t ' 
                          + str(t) + ' cutid ' + str(cutid_limit) + '\n' )
                if from_or_to == 'f':
                    log.joint(' LHS coeff ' + ' Pft ' + str(coeff_P) 
                              + ' Qft ' + str(coeff_Q) + '\n')
                elif from_or_to == 't':
                    log.joint(' LHS coeff ' + ' Ptf ' + str(coeff_P) 
                              + ' Qtf ' + str(coeff_Q) + '\n')
            
            
            if all_data['limit_validity']:
                if from_or_to == 'f':
                    sol_P           = all_data['sol_Pfvalues'][branch]
                    sol_Q           = all_data['sol_Qfvalues'][branch]
                elif from_or_to == 't':
                    sol_P           = all_data['sol_Ptvalues'][branch]
                    sol_Q           = all_data['sol_Qtvalues'][branch]
                sol_violation    = coeff_P * sol_P + coeff_Q * sol_Q - 1

                if sol_violation > FeasibilityTol:
                    log.joint(' WARNING, the limit-envelope cut associated to branch ' 
                              + str(branchid) + ' f ' + str(f) + ' t ' 
                              + str(t) + ' is violated by the AC solution!\n')
                    log.joint(' violation ' + str(sol_violation) + '\n')
                    if from_or_to == 'f':
                        log.joint(' values (AC solution) ' + ' Pft ' 
                                  + str(sol_P) + ' Qft ' + str(sol_Q) + '\n')
                    elif from_or_to == 't':
                        log.joint(' values (AC solution) ' + ' Ptf ' 
                                  + str(sol_P) + ' Qtf ' + str(sol_Q) + '\n')
                    #breakexit('check!')
                else:
                    log.joint(' AC solution satisfies limit cut at branch ' 
                              + str(branchid) + ' with slack ' + str(sol_violation) 
                              + '\n')


            cutexp     = LinExpr()
            if from_or_to == 'f':
                constrname = "limit_cut_"+str(cutid_limit)+"_"+str(branchid)+"r_"+str(rnd)+"_"+str(f)+"_"+str(t)
                cutexp    += coeff_P * Pvar_f[branch] + coeff_Q * Qvar_f[branch]
            elif from_or_to == 't':
                constrname = "limit_cut_"+str(cutid_limit)+"_"+str(branchid)+"r_"+str(rnd)+"_"+str(f)+"_"+str(t)
                cutexp    += coeff_P * Pvar_t[branch] + coeff_Q * Qvar_t[branch]

            themodel.addConstr(cutexp <= 1, name = constrname)
            linenum += 1
    

    all_data['num_jabr_cuts']  = numjabrcuts
    all_data['num_i2_cuts']    = numi2cuts
    all_data['num_limit_cuts'] = numlimitcuts
    all_data['ID_jabr_cuts']   = numjabrcuts
    all_data['ID_i2_cuts']     = numi2cuts
    all_data['ID_limit_cuts']  = numlimitcuts

    all_data['jabr_cuts']       = jabr_cuts
    all_data['jabr_cuts_info']  = jabr_cuts_info
    all_data['i2_cuts']         = i2_cuts
    all_data['i2_cuts_info']    = i2_cuts_info
    all_data['limit_cuts']      = limit_cuts
    all_data['limit_cuts_info'] = limit_cuts_info


def add_i2(log,all_data):
        
    log.joint('\n')
    log.joint(' **** new i2 vars ****\n')

    themodel       = all_data['themodel']
    buses          = all_data['buses']
    branches       = all_data['branches']    
    i2var_f        = all_data['i2var_f']
    newi2          = all_data['newi2']

    for branch in newi2:

      branchcount = branch.count
      f           = branch.f
      t           = branch.t
      count_of_f  = IDtoCountmap[f]
      count_of_t  = IDtoCountmap[t]
      bus_f       = buses[count_of_f]

      upperbound_f = branch.limit**2 / (bus_f.Vmin * bus_f.Vmin)

      i2var_f[branch] = themodel.addVar(obj = 0.0, lb = 0, ub = upperbound_f ,
                                        name = "i2_" + str(branch.count) + "_"
                                        + str(f) + "_" + str(t))

    
    all_data['i2var_f']  = i2var_f
    all_data['themodel'] = themodel

    themodel.update()
    log.joint(' model updated with new i2 variables\n')


