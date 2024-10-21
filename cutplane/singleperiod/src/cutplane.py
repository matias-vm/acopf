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

# For a detailed description of the model variables and constraints please
# check 'cutplane_mtp.py' in '../../multiperiod/src/', which corresponds
# to the multi-period version of this model.

import sys
import math
from log import danoLogger
from gurobipy import *
import numpy as np
import bisect
from myutils import breakexit
import reader
import time
import math
from cuts import *
import random
import os
import platform
import psutil

def gocutplane(log, all_data):

  formulation_start = time.time()
  themodel          = Model("Cutplane")
  buses             = all_data['buses']
  numbuses          = all_data['numbuses']
  branches          = all_data['branches']
  numbranches       = all_data['numbranches']
  gens              = all_data['gens']
  IDtoCountmap      = all_data['IDtoCountmap']

  ############################ LOAD SOLUTION ##################################

  if all_data['ampl_sol']:
    getsol_ampl(log,all_data)

  ################################ VARIABLES ##################################

  cvar    = {}
  svar    = {}
  Pvar_f  = {}
  Qvar_f  = {}
  Pvar_t  = {}
  Qvar_t  = {}
  Pinjvar = {}
  Qinjvar = {}
  GenPvar = {}
  GenQvar = {}
  GenTvar = {}

  log.joint(' creating variables...\n')

  varcount = 0

  #cbus, bus-injection, GenP, GenQ, and GenT variables
  for bus in buses.values():

    maxprod = bus.Vmax*bus.Vmax
    minprod = bus.Vmin*bus.Vmin
          
    ubound = maxprod
    lbound = minprod

    Pubound, Plbound, Qubound, Qlbound = computebalbounds(log, all_data, bus)  
          
    cvar[bus] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                name = "c_" + str(bus.nodeID) + "_" 
                                + str(bus.nodeID))
    Pinjvar[bus] = themodel.addVar(obj = 0.0, lb = Plbound, ub = Pubound, 
                                   name = "IP_"+str(bus.nodeID))
    Qinjvar[bus] = themodel.addVar(obj = 0.0, lb = Qlbound, ub = Qubound, 
                                   name = "IQ_"+str(bus.nodeID))
    
    varcount += 3

    for genid in bus.genidsbycount:
      gen = gens[genid]

      lower = gen.Pmin*gen.status
      upper = gen.Pmax*gen.status

      GenPvar[gen] = themodel.addVar(obj = 0.0, lb = lower, ub = upper, 
                                     name = "GP_" + str(gen.count) + "_" 
                                     + str(gen.nodeID))
      lower = gen.Qmin*gen.status
      upper = gen.Qmax*gen.status

      if bus.nodetype == 3:
        upper = GRB.INFINITY
        lower = -GRB.INFINITY

      GenQvar[gen] = themodel.addVar(obj = 0.0, lb = lower, ub = upper, 
                                     name = "GQ_" + str(gen.count) + "_" 
                                     + str(gen.nodeID))

      varcount += 2



  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]
    maxprod     = buses[count_of_f].Vmax*buses[count_of_t].Vmax
    minprod     = buses[count_of_f].Vmin*buses[count_of_t].Vmin

    ubound      = maxprod
    lbound      = -maxprod
    maxanglerad = branch.maxangle_rad
    minanglerad = branch.minangle_rad

    # Cosine                                                                                                                                           
    if maxanglerad <= 0.5*math.pi:
      # In this case minangle <= 0                                                                                                                   
      if minanglerad >= -0.5*math.pi:
        lbound = minprod*min(math.cos(maxanglerad), math.cos(minanglerad))
      elif minanglerad >= -math.pi:
        lbound = maxprod*math.cos(minangle_rad)  # Which is negative                                                                               
      elif minanglerad >= -1.5*math.pi:
        lbound = -maxprod
      else:
        lbound = -maxprod

    elif maxanglerad <= math.pi:
      if minanglerad >= -0.5*math.pi:
        lbound = maxprod*math.cos(maxanglerad)
      elif minanglerad >= -math.pi:
        lbound = maxprod*min(math.cos(maxanglerad), math.cos(minanglerad))
      elif minanglerad >= -1.5*math.pi:
        lbound = -maxprod
      else:
        lbound = -maxprod

    elif maxanglerad <= 1.5*math.pi:
      lbound = -maxprod    

    elif maxanglerad <= 2*math.pi:
      lbound = -maxprod

    else:
      ubound = maxprod
      lbound = -maxprod

    cvar[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                   name = "c_" + str(branchcount) + "_" 
                                   + str(f) + "_" + str(t))

    #s variables

    if maxanglerad <= math.pi/2:
      ubound = maxprod*math.sin(maxanglerad)

      if  minanglerad >= -0.5*math.pi:
        lbound = maxprod*math.sin(minanglerad)
      elif  minanglerad >= -math.pi:
        lbound = -maxprod
      elif  minanglerad >= -1.5*math.pi:
        ubound = maxprod*max( math.sin(maxanglerad), math.sin(minanglerad))
        lbound = -maxprod
      else:
        ubound = maxprod
        lbound = -maxprod

    elif maxanglerad <= math.pi:
      ubound = maxprod

      if minanglerad >= -0.5*math.pi:
        lbound = maxprod*math.sin(minanglerad)
      elif minanglerad >= -math.pi:
        lbound = -maxprod
      elif minanglerad >= -1.5*math.pi:
        lbound = -maxprod
      else:
        lbound = -maxprod

    elif maxanglerad <= 1.5*math.pi:
      ubound = maxprod

      if minanglerad >= -0.5*math.pi:
        lbound = maxprod*min(math.sin(maxanglerad), math.sin(minanglerad))
      elif minanglerad >= -math.pi:
        lbound = -maxprod
      elif minanglerad >= -1.5*math.pi:
        lbound = -maxprod
      else:
        lbound = -maxprod
    else:
      ubound = maxprod
      lbound = -maxprod

    svar[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                   name = "s_" + str(branchcount) + "_" 
                                   + str(f) + "_" + str(t))

    varcount +=2
    

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]

    ubound = branch.limit 
    lbound = -branch.limit

    Pvar_f[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                     name = "P_" + str(branch.count) + "_" 
                                     + str(f) + "_" + str(t))
    Pvar_t[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                     name = "P_" + str(branch.count) + "_" 
                                     + str(t) + "_" + str(f))
    Qvar_f[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                     name = "Q_" + str(branch.count) + "_" 
                                     + str(f) + "_" + str(t))
    Qvar_t[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                     name = "Q_" + str(branch.count) + "_" 
                                     + str(t) + "_" + str(f))

    varcount +=4 

  #i2 variables
  if all_data['i2']:
    i2var_f = {}

    for branch in branches.values():
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
    
      varcount += 1
      
  themodel.update()

  all_data['themodel']    = themodel
  all_data['cvar']        = cvar
  all_data['svar']        = svar
  all_data['GenPvar']     = GenPvar
  all_data['GenTvar']     = GenTvar
  all_data['Pvar_f']      = Pvar_f
  all_data['Pvar_t']      = Pvar_t
  all_data['Qvar_f']      = Qvar_f
  all_data['Qvar_t']      = Qvar_t

  if all_data['i2']:
    all_data['i2var_f']   = i2var_f
    #all_data['i2var_t']  = i2var_f

  ############################## OBJECTIVE ####################################

  log.joint(' creating the objective...\n')

  # Constant term (fixed costs)                                                                                                                                                            
  constexpr   = LinExpr()
  constobjval = 0

  for gen in gens.values():
    if gen.status > 0:
      constobjval += gen.costvector[gen.costdegree]

  constvar   = themodel.addVar(lb = 1.0, ub = 1.0,name = "constant")
  varcount  += 1
  constexpr += constobjval * constvar

  # Linear and Quadratic terms                                                                                                                                                             
  lincostexpr = LinExpr()
  qcostexpr   = QuadExpr()

  for gen in gens.values():
    lincoeff = gen.costvector[gen.costdegree-1]
    lincostexpr += lincoeff * GenPvar[gen]

    if gen.costdegree == 2 and gen.costvector[0] != 0:
      qcostexpr += gen.costvector[0]*GenPvar[gen]*GenPvar[gen]

  log.joint('   %d variables added\n' %varcount)
      
  themodel.setObjective(constexpr + lincostexpr + qcostexpr)

  log.joint('  objective added\n')
  
  themodel.update()

  ############################# CONSTRAINTS ###################################

  log.joint(' creating constraints...\n')

  constrcount = 0
  count       = 0

  #definition flow variables
  log.joint('  active power flow variables definition\n')

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]

    if branch.status == 0:
      log.joint(' branch ' + str(branch.count) + ' f ' + str(f) + ' t ' 
                + str(t) + ' is OFF\n')
      breakexit('check, reader does not include branches with status = 0')

    #  Gff cff + Gft cft + Bft sft
    constrname = "Pdef_"+str(branch.count)+"_"+str(f)+"_"+str(t)
    expr = LinExpr()
    expr += branch.Gff*cvar[buses[count_of_f]]
    expr += branch.Gft*cvar[branch]
    expr += branch.Bft*svar[branch]
    
    themodel.addConstr(expr == Pvar_f[branch], name = constrname)

    #  Gtt ctt + Gtf cft + Btf stf = Gtt ctt + Gtf cft - Btf sft
    constrname = "Pdef_"+str(branch.count)+"_"+str(t)+"_"+str(f)
    expr = LinExpr()
    expr += branch.Gtt*cvar[buses[count_of_t]]
    expr += branch.Gtf*cvar[branch]
    expr += -branch.Btf*svar[branch]

    themodel.addConstr(expr == Pvar_t[branch], name = constrname)
    
    constrcount += 2
    count       += 2

  log.joint('   %d active power flow definition constraints added\n'%count)

  log.joint('  reactive power flow variables definition\n')
  count = 0

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    constrname = "Qdef_"+str(branch.count)+"_"+str(f)+"_"+str(t)
    
    # -Bff cff - Bft cft + Gft sft
    expr = LinExpr()
    expr += -branch.Bff*cvar[buses[count_of_f]]
    expr += -branch.Bft*cvar[branch]
    expr += +branch.Gft*svar[branch]

    themodel.addConstr(expr == Qvar_f[branch], name = constrname)

    # -Btt ctt - Btf cft + Gtf stf = -Btt ctt - Btf cft - Gtf sft 
    constrname = "Qdef_"+str(branch.count)+"_"+str(t)+"_"+str(f)
    expr = LinExpr()
    expr += -branch.Btt*cvar[buses[count_of_t]]
    expr += -branch.Btf*cvar[branch]
    expr += -branch.Gtf*svar[branch]

    themodel.addConstr(expr == Qvar_t[branch], name = constrname)

    constrcount += 2
    count       += 2
  log.joint('   %d reactive power flow definition constraints added\n'%count)

  #balance constraints
  log.joint('  active power injection constraints\n')
  count = 0

  for bus in buses.values():
    constrname = "PBaldef"+str(bus.nodeID)
    expr = LinExpr()
    for branchid in bus.frombranchids.values():
      expr += Pvar_f[ branches[branchid] ]
    for branchid in bus.tobranchids.values():
      expr += Pvar_t[ branches[branchid] ]

    if ( (bus.Gs != 0) and ( ( len(bus.frombranchids) != 0 ) 
                             or ( len(bus.tobranchids) != 0 ) ) ):
      expr += bus.Gs*cvar[bus]

    themodel.addConstr(expr == Pinjvar[bus], name = constrname)

    constrcount += 1
    count       += 1
  
  log.joint('   %d active power injection constraints added\n'%count)
  log.joint('  reactive power injection constraints\n')
  count = 0

  for bus in buses.values():
    constrname = "QBaldef"+str(bus.nodeID)
    expr = LinExpr()
    for branchid in bus.frombranchids.values():
      expr += Qvar_f[ branches[branchid] ]
    for branchid in bus.tobranchids.values():
      expr += Qvar_t[ branches[branchid] ]
 
    if ( (bus.Bs != 0) and ( ( len(bus.frombranchids) != 0 ) 
                             or ( len(bus.tobranchids) != 0 ) ) ):
      expr += (-bus.Bs)*cvar[bus]

    themodel.addConstr(expr == Qinjvar[bus], name = constrname)
    
    constrcount += 1
    count       += 1
  log.joint('   %d reactive power injection constraints added\n'%count)

  #definition bus-injection variables

  log.joint('  adding injection definition constraints...\n')

  count = 0

  for bus in buses.values():
    constrname = "Bus_PInj_"+str(bus.nodeID)
    expr = LinExpr()

    if len(bus.genidsbycount) > 0:
      for genid in bus.genidsbycount:
        gen = gens[genid]
        expr += GenPvar[gen]

    themodel.addConstr(Pinjvar[bus] == expr - bus.Pd, name = constrname)

    constrname = "Bus_QInj_"+str(bus.nodeID)
    expr = LinExpr()

    if len(bus.genidsbycount) > 0:
      for genid in bus.genidsbycount:
        gen = gens[genid]
        expr += GenQvar[gen]

    themodel.addConstr(Qinjvar[bus] == expr - bus.Qd, name = constrname)

    constrcount += 2
    count       += 2

  log.joint('   %d power injection definitions added\n'%count)
  
  #definition i2 variables
  if all_data['i2']:
    constrcount += i2_def(log,all_data)

  #jabr inequalities
  if all_data['jabr_inequalities']:
    constrcount += jabr_inequalities(log,all_data)

  #i2 inequalities
  if all_data['i2_inequalities']:
    constrcount += i2_inequalities(log,all_data)

  #limit constraints
  if all_data['limit_inequalities']:
    constrcount += limit_inequalities(log,all_data)
  
  log.joint('  %d constraints added\n'%constrcount)
    
  themodel.update()

  formulation_end = time.time()

  all_data['formulation_time'] = formulation_end - formulation_start

  log.joint(' formulation time: %g\n' % all_data['formulation_time'])

  # Write model to .lp file
  #if all_data['writelps']:
  #  log.joint(' writing to lpfile ' + all_data['lpfilename'] + '\n')  
  #  themodel.write(all_data['lpfilename'])

    #breakexit('writelp without cuts, then remove break')
    
  ###################### INIT DATA STRUCTURES FOR CUTS ########################

  if all_data['jabrcuts']:
    jabr_cuts_info = all_data['jabr_cuts_info']

    for branch in branches.values():
      jabr_cuts_info[branch] = {}

  if all_data['i2cuts']:
    i2_cuts_info = all_data['i2_cuts_info']

    for branch in branches.values():
      i2_cuts_info[branch] = {}

  if all_data['limitcuts']:
    limit_cuts_info = all_data['limit_cuts_info']

    for branch in branches.values():
      limit_cuts_info[branch] = {}

  ######################## FIXING/WRITING AN AC SOLUTION ######################

  if all_data['fixflows']:
    fixflows(log,all_data)
    if all_data['fixcs'] == 0:
      return None

  if all_data['fixcs']:
    fixcs(log,all_data)
    return None

  ########################### SOLVER PARAMETERS ###############################

  themodel.Params.Method    = all_data['solver_method']
  themodel.Params.Crossover = all_data['crossover'] 

  if all_data['solver_method'] == 2:
    themodel.Params.BarHomogeneous = 1
    themodel.Params.BarConvTol     = 1e-6 # Hard coded tolerances
    themodel.Params.FeasibilityTol = 1e-6
    themodel.Params.OptimalityTol  = 1e-6

  themodel.Params.NumericFocus = 1
  themodel.Params.OutPutFlag = 1
  
  ######################### READING AND LOADING CUTS ##########################

  if all_data['addcuts']:

    t0_cuts = time.time()

    if all_data['max_rounds'] > 1:
      add_cuts_ws(log,all_data)
    else:
      add_cuts(log,all_data)

    themodel.update()

    t1_cuts = time.time()

    all_data['addcuts_time'] = t1_cuts - t0_cuts

    log.joint(' pre-computed cuts added and model updated\n')

    log.joint(' reading and loading cuts time = '
              + str(all_data['addcuts_time']) + '\n')

    if all_data['writelps']:
      themodel.write(all_data['casename']+'_precomputed_cuts.lp')
      log.joint(' model with precomputed written to .lp file\n')
      #breakexit('precomputed cuts')
      
  ########################## CUTPLANE MAIN LOOP ###############################

  all_data['round']                  = 1
  all_data['runtime']                = time.time() - all_data['T0']
  all_data['cumulative_solver_time'] = 0
  all_data['ftol_counter']           = 0
  oldobj                             = 1
  gap                                = 1e20

  while ((all_data['round'] <= all_data['max_rounds']) and 
         (all_data['runtime'] <= all_data['max_time']) and 
         (all_data['ftol_counter'] <= all_data['ftol_iterates'])):
    
  
    ############################ SOLVING MODEL ################################

    cutplane_optimize(log,all_data)

    ########################### STORING SOLUTION ##############################

    log.joint(' storing current solution ...\n')

    all_data['Pfvalues']   = themodel.getAttr("X",Pvar_f)
    all_data['Qfvalues']   = themodel.getAttr("X",Qvar_f)
    all_data['Ptvalues']   = themodel.getAttr("X",Pvar_t)
    all_data['Qtvalues']   = themodel.getAttr("X",Qvar_t)
    all_data['cvalues']    = themodel.getAttr("X",cvar)
    all_data['svalues']    = themodel.getAttr("X",svar)
    all_data['GenPvalues'] = themodel.getAttr("X",GenPvar)
    all_data['GenQvalues'] = themodel.getAttr("X",GenQvar)

    if all_data['i2']:
      all_data['i2fvalues'] = themodel.getAttr("X",i2var_f)

    log.joint(' done storing values\n')
     
    ########################## CHECK OBJ IMPROVEMENT ##########################

    if ((all_data['objval'] - oldobj)/oldobj) < all_data['ftol']:
      all_data['ftol_counter'] += 1
    else:
      all_data['ftol_counter'] = 0

    oldobj              = all_data['objval']
    all_data['runtime'] = time.time() - all_data['T0']

    ########################### ROUND STATISTICS ##############################

    cutplane_stats(log,all_data)

    ######################### SUMMARY EXPERIMENTS #############################

    avg = 0
    numbuses = 0

    for bus in buses.values():
      nodeID            = bus.nodeID
      constrname        = "PBaldef"+str(bus.nodeID)
      constr            = themodel.getConstrByName(constrname)
      dual              = - constr.Pi/all_data['baseMVA']
      avg              += dual
      numbuses += 1
      log.joint(constrname + ' = ' + str(dual) + '\n')
    
    log.joint(' welfare ' + str(all_data['objval']) + '\n')
    log.joint(' avg price ' + str(avg/numbuses) + '\n')
    breakexit('check prices')

    ###########################################################################
    all_data['runtime'] = time.time() - all_data['T0']

    log.joint("\n writing casename, opt stauts, obj and runtime to summary_ws.log\n")

    summary_ws = open("summary_ws.log","a+") 

    numcutsadded = 0
    numcuts      = 0
    if all_data['jabrcuts']:
      numcutsadded += all_data['ID_jabr_cuts']
      numcuts      += all_data['num_jabr_cuts']
    if all_data['i2cuts']:
      numcutsadded += all_data['ID_i2_cuts']
      numcuts      += all_data['num_i2_cuts']
    if all_data['limitcuts']:
      numcutsadded += all_data['ID_limit_cuts']
      numcuts      += all_data['num_limit_cuts']
        
    summary_ws.write(' case ' + all_data['casename'] + ' opt_status ' 
                     + str(all_data['optstatus']) + ' obj ' 
                     + str(all_data['objval']) + ' runtime ' 
                     + str(all_data['runtime']) + ' iterations ' 
                     + str(all_data['round']) + ' rndcuts '
                     + str((all_data['round']-1)) + ' numcutsadded '
                     + str(numcutsadded) + ' numcuts '
                     + str(numcuts) +  '\n')

    summary_ws.close()
    
    ############################ TERMINATION #################################

    if (all_data['round'] >= all_data['max_rounds']):

      if all_data['writesol']:
        writesol(log,all_data)
        writesol_allvars(log,all_data)
      
      summary_ws = open("summary_ws.log","a+") 
      summary_ws.write(' rounds limit reached!\n\n')
      summary_ws.close()
      log.joint(' rounds limit reached!\n')
      log.joint(' bye\n')
      return None
          
    if all_data['runtime'] >= all_data['max_time']:

      if all_data['writesol']:
        writesol(log,all_data)
        writesol_allvars(log,all_data)
        
      summary_ws = open("summary_ws.log","a+") 
      summary_ws.write(' time limit reached!\n\n')
      summary_ws.close()
      log.joint(' time limit reached!\n')
      log.joint(' bye\n')
      return None

    if (all_data['ftol_counter'] >= all_data['ftol_iterates']):

      if all_data['writesol']:
        writesol(log,all_data)
        writesol_allvars(log,all_data)
      
      summary_ws = open("summary_ws.log","a+") 
      summary_ws.write(' poor consecutive obj improvement limit reached!\n\n')
      summary_ws.close()
      log.joint(' poor consecutive obj improvement limit reached\n')
      log.joint(' bye\n')
      return None

    ##########################################################################

    if all_data['max_rounds'] == 1:

      if all_data['writesol']:
        writesol(log,all_data)      
        writesol_allvars(log,all_data)
        
      log.joint(' bye!\n')
      return None

    ############################### CUTS ######################################

    # Cut computations and management
    cutplane_cuts(log,all_data)

    # Cut statistics
    cutplane_cutstats(log,all_data)
    
    themodel.update()

    log.joint(' model updated\n')
    log.joint('\n')

    ############################### WRITE CUTS ################################

    if all_data['writecuts']:
      write_cuts(log,all_data)

    ############################### WRITE LPS #################################
    

    if all_data['writelps'] and ( all_data['round'] > 0 ):
      name = 'post_cuts' + '_' + str(all_data['round']) + '.lp'
      themodel.write(name)
      log.joint(' model with new cuts written to .lp file\n')


    all_data['round'] += 1


    ###########################################################################

def fixflows(log,all_data):

  if (all_data['ampl_sol'] == 0):
    log.joint(' cannot fix flows since no solution has been loaded!\n')
    return None

  themodel    = all_data['themodel']
  branches    = all_data['branches']
  tolerance   = all_data['tol_fix']
  Pvar_f      = all_data['Pvar_f']
  Pvar_t      = all_data['Pvar_t']
  Qvar_f      = all_data['Qvar_f']
  Qvar_t      = all_data['Qvar_t']

  sol_Pfvalues = all_data['sol_Pfvalues']
  sol_Ptvalues = all_data['sol_Ptvalues']
  sol_Qfvalues = all_data['sol_Qfvalues']
  sol_Qtvalues = all_data['sol_Qtvalues']

  for branch in branches.values():
        
    sol_Pf = sol_Pfvalues[branch]
    sol_Pt = sol_Ptvalues[branch]
    sol_Qf = sol_Qfvalues[branch]
    sol_Qt = sol_Qtvalues[branch]

    #Pf
    ubound_Pf = sol_Pf + tolerance
    lbound_Pf = sol_Pf - tolerance

    Pvar_f[branch].setAttr("ub",ubound_Pf)
    Pvar_f[branch].setAttr("lb",lbound_Pf)

    #Pt
    ubound_Pt = sol_Pt + tolerance
    lbound_Pt = sol_Pt - tolerance

    Pvar_t[branch].setAttr("ub",ubound_Pt)
    Pvar_t[branch].setAttr("lb",lbound_Pt)

    #Qf
    ubound_Qf = sol_Qf + tolerance
    lbound_Qf = sol_Qf - tolerance

    Qvar_f[branch].setAttr("ub",ubound_Qf)
    Qvar_f[branch].setAttr("lb",lbound_Qf)

    #Qt
    ubound_Qt = sol_Qt + tolerance
    lbound_Qt = sol_Qt - tolerance

    Qvar_t[branch].setAttr("ub",ubound_Qt)
    Qvar_t[branch].setAttr("lb",lbound_Qt)

  themodel.update()
  themodel.write('fixflows.lp')
  log.joint('check fixflows.lp\n')  


def fixcs(log,all_data):

  if (all_data['ampl_sol'] == 0):
    log.joint(' cannot fix flows since no solution has been loaded!\n')
    return None

  themodel   = all_data['themodel']
  tolerance  = all_data['tol_fix']
  buses      = all_data['buses']
  branches   = all_data['branches']
  cvar       = all_data['cvar']
  svar       = all_data['svar']
  sol_cvalues = all_data['sol_cvalues']
  sol_svalues = all_data['sol_svalues']

  for bus in buses.values():

    sol_v2 = all_data['sol_cvalues'][bus]

    ubound_v2 = sol_v2 + tolerance
    lbound_v2 = sol_v2 - tolerance

    cvar[bus].setAttr("ub",ubound_v2)
    cvar[bus].setAttr("lb",lbound_v2)

  for branch in branches.values():

    sol_c = all_data['sol_cvalues'][branch]
    sol_s = all_data['sol_svalues'][branch]

    #c
    ubound_c = sol_c + tolerance
    lbound_c = sol_c - tolerance

    cvar[branch].setAttr("ub",ubound_c)
    cvar[branch].setAttr("lb",lbound_c)

    #s
    ubound_s = sol_s + tolerance
    lbound_s = sol_s - tolerance

    svar[branch].setAttr("ub",ubound_s)
    svar[branch].setAttr("lb",lbound_s)

  themodel.update()
  themodel.write('fixCS.lp')
  log.joint('check fixCS.lp\n')


def computebalbounds(log, all_data, bus):

  # We get max/min generations
  
  baseMVA = all_data['baseMVA']
  gens = all_data['gens']

  Pubound = Plbound = 0
  Qubound = Qlbound = 0

  for gencounter in bus.genidsbycount:
    if gens[gencounter].status:
      Pubound += gens[gencounter].Pmax
      Plbound += gens[gencounter].Pmin
      Qubound += gens[gencounter].Qmax
      Qlbound += gens[gencounter].Qmin

      if bus.nodetype == 3:
        log.joint('IQ bound for generator %d buscnt %d (reference) set to infinities.\n'%(gencounter, bus.count))
        Qubound = GRB.INFINITY
        Qlbound = -GRB.INFINITY        

  
  Pubound -= bus.Pd
  Plbound -= bus.Pd
  Qubound -= bus.Qd
  Qlbound -= bus.Qd

  if bus.nodetype == 4:
    Pubound = Plbound = Qubound = Qlbound = 0

  
  return Pubound, Plbound, Qubound, Qlbound


def computeangles(log,all_data):
  
  cvalues  = all_data['cvalues']
  svalues  = all_data['svalues']
  branches = all_data['branches']

  for branch in branches.values():
    s = svalues[branch]
    c = cvalues[branch]

    log.joint(' branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' angle ' + str(math.atan(s/c)) + '\n')

  

def computei2value(log,all_data,branch,mp_c,mp_s,mp_cbusf,mp_cbust):

  ratio  = branch.ratio
  y      = branch.y
  g      = y.real
  b      = y.imag
  bshunt = branch.bc
  angle  = branch.angle_rad

  i2f  = 0
  i2f += (g*g + b*b)/(ratio*ratio) * ( (mp_cbusf/(ratio*ratio)) + mp_cbust - (2/ratio) * ( mp_c * math.cos(angle) + mp_s * math.sin(angle) ) )
  i2f += b*bshunt/(ratio**3) * ( (mp_cbusf/ratio) - (mp_c * math.cos(angle) + mp_s * math.sin(angle) ))
  i2f += g*bshunt/(ratio**3) * ( mp_s * math.cos(angle) - mp_c * math.sin(angle) )
  i2f += (bshunt*bshunt*mp_cbusf/(4*(ratio**4)) )

  return i2f
  

def jabr_inequalities(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  buses          = all_data['buses']
  cvar           = all_data['cvar']
  svar           = all_data['svar']
  IDtoCountmap   = all_data['IDtoCountmap']
  FeasibilityTol = all_data['tolerance']

  if all_data['ampl_sol'] and all_data['jabr_validity']:
    log.joint('  adding and checking validity of Jabr inequalities wrt a solution\n')
    log.joint('   feasibility tolerance ' + str(FeasibilityTol) + '\n')
  else:
    log.joint('  Jabr inequalities\n')

  maxviolation = 0
  violated     = 0
  maxbranch    = -1
  maxbusf      = -1
  maxbust      = -1
  counter_jabr = 0

  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]

    if all_data['ampl_sol'] and all_data['jabr_validity']:
      sol_c         = all_data['sol_cvalues'][branch]
      sol_s         = all_data['sol_svalues'][branch]
      sol_cbusf     = all_data['sol_cvalues'][buses[count_of_f]]
      sol_cbust     = all_data['sol_cvalues'][buses[count_of_t]]
      relviolation = violation = sol_c * sol_c + sol_s * sol_s - sol_cbusf * sol_cbust
      #relviolation = violation / ( ( coeff_Pft**2 + coeff_Qft**2 + coeff_cff**2 + coeff_i2ft**2 )**0.5 )
      if relviolation > maxviolation:
        maxviolation = relviolation
        maxbranch    = branch.count
        maxbusf      = f
        maxbust      = t

      if relviolation > FeasibilityTol:
        violated += 1
        log.joint('   WARNING, the Jabr inequality associated to branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' is violated by the AC solution!\n')
        log.joint('   violation ' + str(violation) + '\n')
        log.joint('   relative violation ' + str(relviolation) + '\n')
        log.joint('   values (AC solution) ' + ' cft ' + str(sol_c) + ' sft ' + str(sol_s) + ' cff ' + str(sol_cbusf) + ' ctt ' + str(sol_cbust) + '\n' )
        #breakexit('check!')
      else:
        log.joint('   AC solution satisfies Jabr inequality at branch ' + str(branchcount) + ' with slack ' + str(relviolation) + '\n')

    counter_jabr += 1
    trigexp       = QuadExpr()
    constrname    = "jabr_"+str(branchcount)+"_"+str(f)+"_"+str(t)
    trigexp      += cvar[branch]*cvar[branch] + svar[branch]*svar[branch] - cvar[buses[count_of_f]]*cvar[buses[count_of_t]]

    themodel.addConstr(trigexp <= 0, name = constrname)

  if all_data['ampl_sol'] and all_data['jabr_validity']:
    log.joint('  max violation of Jabr-inequalities by AC solution ' + str(maxviolation) + ' at branch ' + str(maxbranch) + ' f ' + str(maxbusf) + ' t ' + str(maxbust) + '\n')
    log.joint('  number of violated Jabr-inequalities ' + str(violated) + '\n')
    #breakexit('  check Jabr violation')

  log.joint('   %d Jabr inequalities added\n'%counter_jabr) 

  return counter_jabr



def i2_def(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  buses          = all_data['buses']
  i2var_f        = all_data['i2var_f']
  cvar           = all_data['cvar']
  svar           = all_data['svar']
  IDtoCountmap   = all_data['IDtoCountmap']
  
  counter_i2def  = 0

  log.joint('  i2 variables definition\n')

  for branch in branches.values():
    expr_f = LinExpr()
    #expr_t = LinExpr()

    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]
    ratio       = branch.ratio
    y           = branch.y
    g           = y.real
    b           = y.imag
    bshunt      = branch.bc
    angle       = branch.angle_rad

    expr_f += (g*g + b*b)/(ratio*ratio) * ( (cvar[buses[count_of_f]]/(ratio*ratio)) + cvar[buses[count_of_t]] - (2/ratio) * ( cvar[branch] * math.cos(angle) + svar[branch] * math.sin(angle) ) )
    expr_f += b*bshunt/(ratio**3) * ( (cvar[buses[count_of_f]]/ratio) - (cvar[branch] * math.cos(angle) + svar[branch] * math.sin(angle) )) #since gsh is 0
    expr_f += g*bshunt/(ratio**3) * ( svar[branch] * math.cos(angle) - cvar[branch] * math.sin(angle) )
    expr_f += (bshunt*bshunt*cvar[buses[count_of_f]]/(4*(ratio**4)) )
    
    counter_i2def  += 1

    constrname_f    = 'i2def_'+str(branch.count)+"_"+str(f) + "_" + str(t)
    themodel.addConstr(expr_f == i2var_f[branch],name = constrname_f) 


  log.joint('   %d i2 definition constraints added\n'%counter_i2def) 

  return counter_i2def


def i2_inequalities(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  buses          = all_data['buses']
  cvar           = all_data['cvar']
  i2var_f        = all_data['i2var_f']
  Pvar_f         = all_data['Pvar_f']
  Qvar_f         = all_data['Qvar_f']
  IDtoCountmap   = all_data['IDtoCountmap']
  FeasibilityTol = all_data['tolerance']

  if all_data['ampl_sol'] and all_data['i2_validity']:
    log.joint('  adding and checking validity of i2 inequalities wrt a solution\n')
    log.joint('   feasibility tolerance ' + str(FeasibilityTol) + '\n')
  else:
    log.joint('  i2 inequalities\n')

  maxviolation = 0
  violated     = 0
  maxbranch    = -1
  maxbusf      = -1
  maxbust      = -1
  maxi2f       = -1
  maxcff       = -1
  maxPf        = -1
  maxQf        = -1

  counter_i2   = 0

  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]

    if all_data['ampl_sol'] and all_data['i2_validity']:
      sol_Pf        = all_data['sol_Pfvalues'][branch]
      sol_Qf        = all_data['sol_Qfvalues'][branch]
      sol_c         = all_data['sol_cvalues'][branch]
      sol_s         = all_data['sol_svalues'][branch]
      sol_cbusf     = all_data['sol_cvalues'][buses[count_of_f]]
      sol_cbust     = all_data['sol_cvalues'][buses[count_of_t]]
      sol_i2f       = computei2value(log,all_data,branch,sol_c,sol_s,sol_cbusf,sol_cbust)
      relviolation = violation = sol_Pf * sol_Pf + sol_Qf * sol_Qf - sol_cbusf * sol_i2f
      #relviolation = violation / ( ( coeff_Pft**2 + coeff_Qft**2 + coeff_cff**2 + coeff_i2ft**2 )**0.5 )
      if relviolation > maxviolation:
        maxviolation = relviolation
        maxbranch    = branch.count
        maxbusf      = f
        maxbust      = t
        maxi2f       = sol_i2f
        maxcff       = sol_cbusf
        maxPf        = sol_Pf
        maxQf        = sol_Qf

      if relviolation > FeasibilityTol:
        violated += 1
        log.joint('   WARNING, the i2 inequality associated to branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' is violated by the AC solution!\n')
        log.joint('   violation ' + str(violation) + '\n')
        log.joint('   relative violation ' + str(relviolation) + '\n')
        log.joint('   values (AC solution) ' + ' Pft ' + str(sol_Pf) + ' Qft ' + str(sol_Qf) + ' cff ' + str(sol_cbusf) + ' i2ft ' + str(sol_i2f) + '\n' )
        #breakexit('check!')
      else:
        log.joint('   AC solution satisfies i2 inequality at branch ' + str(branchcount) + ' with slack ' + str(relviolation) + '\n')

    counter_i2 += 1
    trigexp     = QuadExpr()
    constrname  = "i2_"+str(branchcount)+"_"+str(f)+"_"+str(t)
    trigexp    += Pvar_f[branch]**2 + Qvar_f[branch]**2 - cvar[buses[count_of_f]] * i2var_f[branch]
    themodel.addConstr(trigexp <= 0, name = constrname)

  if all_data['ampl_sol'] and all_data['i2_validity']:
    log.joint('  max violation of i2 inequalities by AC solution ' + str(maxviolation) + ' at branch ' + str(maxbranch) + ' f ' + str(maxbusf) + ' t ' + str(maxbust) + '\n')
    log.joint('  values (AC solution) ' + ' Pft ' + str(maxPf) + ' Qft ' + str(maxQf) + ' cff ' + str(maxcff) + ' i2ft ' + str(maxi2f) + '\n' )
    log.joint('  number of violated i2 inequalities ' + str(violated) + '\n')
    #breakexit('  check i2 violation')

  log.joint('   %d i2 inequalities added\n'%counter_i2) 

  return counter_i2  


def limit_inequalities(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  Pvar_f         = all_data['Pvar_f']
  Pvar_t         = all_data['Pvar_t']
  Qvar_f         = all_data['Qvar_f']
  Qvar_t         = all_data['Qvar_t']

  counter_limit  = 0

  log.joint('  limit inequalities\n')

  for branch in branches.values():
    if True:
    #if branch.constrainedflow:
      branchcount = branch.count
      f           = branch.f
      t           = branch.t

      constrname = "limit_f_"+str(branchcount)+"_"+str(f)+"_"+str(t)
      limexp = QuadExpr()
      limexp += Pvar_f[branch]*Pvar_f[branch] + Qvar_f[branch]*Qvar_f[branch]
      themodel.addConstr(limexp <= branch.limit**2, name = constrname)

      constrname = "limit_t_"+str(branchcount)+"_"+str(t)+"_"+str(f)
      limexp = QuadExpr()
      limexp += Pvar_t[branch]*Pvar_t[branch] + Qvar_t[branch]*Qvar_t[branch]
      themodel.addConstr(limexp <= branch.limit**2, name = constrname)

      counter_limit += 2

  log.joint('   %d limit inequalities added\n'%counter_limit) 

  return counter_limit


def cutplane_stats(log,all_data):

  themodel = all_data['themodel']
  
  log.joint('\n ******************** round statistics **********************\n') 
    
  log.joint(' casename = %s\n' % all_data['casename'] )
  log.joint(' round = %g\n' % all_data['round'] )
  log.joint(' objective = %g\n' % all_data['objval'])
  log.joint(' solver status = ' + str(all_data['optstatus']) 
            + ' solver method = ' + str(themodel.params.method) + '\n')
  if themodel.params.method == 2:
    log.joint(' crossover ' + str(themodel.params.crossover) + '\n')
  log.joint(' BarConTol = ' + str(themodel.params.BarConvTol) 
            + ' FeasTol = ' + str(themodel.params.FeasibilityTol) 
            + ' OptTol = ' + str(themodel.params.OptimalityTol) + '\n') 


  if all_data['addcuts'] and all_data['round'] == 1:
    log.joint(' -- precomputed cuts --\n')
    log.joint(' Jabr-envelope cuts = %d\n' 
               % all_data['addcuts_numjabrcuts'])
    log.joint(' i2-envelope cuts = %d\n' 
               % all_data['addcuts_numi2cuts'])
    log.joint(' Limit-envelope cuts = %d\n' 
               % all_data['addcuts_numlimitcuts'])
  else:
    log.joint(' -- cut parameters --\n')
    log.joint(' cut age limit = ' 
              + str(all_data['cut_age_limit']) + '\n')
    log.joint(' parallel-cuts threshold = ' 
              + str(all_data['threshold_dotprod']) + '\n')
    log.joint(' threshold = ' 
              + str(all_data['tolerance']) + '\n') 

    if all_data['jabrcuts']:
      log.joint(' -- Jabr-envelope cuts --\n')
      log.joint(' cuts = %d\n' 
                % all_data['num_jabr_cuts'])
      log.joint(' top percent of most violated cuts added = %g\n' 
                % (100*all_data['most_violated_fraction_jabr']) )
      log.joint(' max error (abs) = %g\n' 
                % all_data['max_error_jabr'])
      log.joint(' cut threshold = %g\n' 
                % all_data['threshold']) 

    if all_data['i2cuts']:
      log.joint(' -- i2-envelope cuts --\n')
      log.joint(' cuts = %d\n' 
                % all_data['num_i2_cuts'])
      log.joint(' top percent of most violated cuts added = %g\n' 
                % (100*all_data['most_violated_fraction_i2']))
      log.joint(' max error (abs) = %g\n'
                % all_data['max_error_i2'])
      log.joint(' cut threshold = %g\n' 
                % all_data['threshold_i2']) 

    if all_data['limitcuts']:
      log.joint(' -- Limit-envelope cuts --\n')
      log.joint(' cuts = %d\n' 
                % all_data['num_limit_cuts'])
      log.joint(' top percent of most violated added = %g\n' 
                % (100*all_data['most_violated_fraction_limit']))
      log.joint(' max error (abs) = %g\n' 
                % all_data['max_error_limit'])
      log.joint(' cut threshold = %g\n' 
                % all_data['threshold_limit'])
    
  log.joint(' -- runtimes --\n')
  log.joint(' solver runtime (current round) = %g\n'
            % all_data['solvertime'])
  log.joint(' cumulative solver time = %g\n' 
            % all_data['cumulative_solver_time'])
  
  timenow = time.time()

  log.joint(' time so far (overall - formulation time) = %g\n'
            % (timenow - all_data['T0'] - all_data['formulation_time']))
  log.joint(' time so far (overall) = %g\n'
            % (timenow - all_data['T0']))
  log.joint(' max runtime = %g\n'
            % all_data['max_time'])
  if all_data['ftol_counter']:
    log.joint(' minimum obj improvement (ftol) = %g\n'
              % all_data['ftol'])
    log.joint(' consecutive rounds with poor obj improvement = %d\n' 
              % all_data['ftol_counter'])

  log.joint(' ************************************************************\n') 



def cutplane_cutstats(log,all_data):

  themodel = all_data['themodel']

  log.joint('\n *********************** cuts statistics ********************\n') 
  if all_data['jabrcuts']:
    log.joint(' -- Jabr-envelope cuts --\n')
    log.joint(' cuts = %d\n'
              % all_data['num_jabr_cuts'])
    if all_data['addcuts']:
      log.joint('  precomputed cuts = %d\n' 
              % all_data['addcuts_numjabrcuts'])
    log.joint('  added in current round = %d\n' 
              % all_data['num_jabr_cuts_added'])
    if all_data['dropjabr']:
      log.joint('  dropped in current round = %d\n' 
                % all_data['num_jabr_cuts_dropped'])
    log.joint(' added (overall) = %d\n'
              % all_data['ID_jabr_cuts'])
    if all_data['dropjabr']:
      log.joint(' dropped (overall) = %d\n'
                % all_data['total_jabr_dropped'])

  if all_data['i2cuts']:
    log.joint(' -- i2-envelope cuts --\n')
    log.joint(' cuts = %d\n'
              % all_data['num_i2_cuts'])
    if all_data['addcuts']:
      log.joint('  precomputed cuts = %d\n' 
                % all_data['addcuts_numi2cuts'])
    log.joint('  added in current round = %d\n' 
              % all_data['num_i2_cuts_added'])
    if all_data['dropi2']:
      log.joint('  dropped in current round = %d\n'
                % all_data['num_i2_cuts_dropped'] )
  
    log.joint(' added (overall) = %d\n'
              % all_data['ID_i2_cuts'])

    if all_data['dropi2']:
      log.joint(' dropped (overall) = %d\n'
                % all_data['total_i2_dropped'])

  if all_data['limitcuts']:
    log.joint(' -- Limit-envelope cuts --\n')
    log.joint(' cuts = %d\n'
              % all_data['num_limit_cuts'])
    if all_data['addcuts']:
      log.joint('  precomputed cuts = %d\n' 
                % all_data['addcuts_numlimitcuts'])
    log.joint('  added in current round = %d\n'
              % all_data['num_limit_cuts_added'])
    if all_data['droplimit']:
      log.joint('  dropped in current round = %d\n'
                % all_data['num_limit_cuts_dropped'] )
    log.joint(' added (overall) = %d\n'
               % all_data['ID_limit_cuts'])
    if all_data['droplimit']:
      log.joint(' dropped (overall) = %d\n'
                % all_data['total_limit_dropped'])

  log.joint(' ************************************************************\n\n') 


def cutplane_cuts(log,all_data):

  log.joint('\n starting cut procedure ...\n')

  t0_cuts = time.time()

  t0_jabr = time.time()

  if all_data['jabrcuts']:
    jabr_cuts(log,all_data)

    if all_data['NO_jabrs_violated']:
      if all_data['threshold'] > all_data['tolerance']:
        all_data['threshold'] *= 1e-01
        log.joint(' threshold updated to ' + str(all_data['threshold']) + '\n' )

  t1_jabr = time.time()

  log.joint(' time spent on Jabr-cuts ' + str(t1_jabr - t0_jabr) + '\n')


  t0_i2 = time.time()

  if all_data['i2cuts']:
    i2_cuts(log,all_data)

    if all_data['NO_i2_cuts_violated']:
      if all_data['threshold_i2'] > all_data['tolerance']:
        all_data['threshold_i2'] *= 1e-01
        log.joint(' threshold updated to ' + str(all_data['threshold_i2']) + '\n' )

  t1_i2 = time.time()
  log.joint(' time spent on i2-cuts ' + str(t1_i2 - t0_i2) + '\n')

  t0_lim = time.time()

  if all_data['limitcuts']:
    limit_cuts(log,all_data)
    if all_data['NO_limit_cuts_violated']:
      if all_data['threshold_limit'] > all_data['tolerance']:
        all_data['threshold_limit'] *= 1e-01
        log.joint(' threshold updated to ' + str(all_data['threshold_limit']) + '\n' )

  t1_lim = time.time()
  log.joint(' time spent on lim-cuts ' + str(t1_lim - t0_lim) + '\n')

  t1_cuts = time.time()

  log.joint('\n time spent on cuts ' + str(t1_cuts - t0_cuts) + '\n')



def cutplane_cutmanagement(log,all_data):

  t0_cutmanag = time.time()

  if all_data['jabrcuts'] and all_data['dropjabr'] and (all_data['round'] >= all_data['cut_age_limit']):
    drop_jabr(log,all_data)

  if all_data['i2cuts'] and all_data['dropi2'] and (all_data['round'] >= all_data['cut_age_limit']):
    drop_i2(log,all_data)

  if all_data['limitcuts'] and all_data['droplimit'] and (all_data['round'] >= all_data['cut_age_limit']):
    drop_limit(log,all_data)

  log.joint('\n')

  t1_cutmanag = time.time()


def cutplane_optimize(log,all_data):

  themodel = all_data['themodel']

  if all_data['jabr_inequalities'] or all_data['limit_inequalities']:
    themodel.params.QCPDual = 1

  log.joint(' solving model with method ' + str(themodel.params.method) + '\n')
  log.joint(' crossover ' + str(themodel.params.crossover) + '\n')

  t0_solve = time.time()
  themodel.optimize()
  t1_solve = time.time()

  if themodel.status == GRB.status.INF_OR_UNBD:
    log.joint(' -> LP infeasible or unbounded\n')
    log.joint(' turning presolve off and reoptimizing\n')
    themodel.params.presolve = 0

    t0_solve = time.time()
    themodel.optimize()
    t1_solve = time.time()

    all_data['runtime'] = t1_solve - all_data['T0']

    log.joint(' writing casename, opt status, and runtime to summary_ws.log\n')

    summary_ws = open("summary_ws.log","a+") 
    summary_ws.write(' case ' + all_data['casename'] + ' opt_status ' 
                     + str(themodel.status) + ' runtime ' 
                     + str(all_data['runtime']) + ' iterations ' 
                     + str(all_data['round']) + '\n')

    summary_ws.close()

    log.joint(' optimization status ' + str(themodel.status) + '\n')
    log.joint(' solver runtime current round = %g\n' % (t1_solve - t0_solve) )
    log.joint(' overall time = %g\n' % all_data['runtime'] )
    log.joint(' bye.\n')
    return None

  elif themodel.status == GRB.status.INFEASIBLE:
    log.joint(' -> LP infeasible\n')

    all_data['runtime'] = time.time() - all_data['T0']

    themodel.computeIIS()
    model.write('model.ilp')

    log.joint(' writing casename, opt status, and runtime to summary_ws.log\n')

    summary_ws = open("summary_ws.log","a+")
    summary_ws.write(' case ' + all_data['casename'] + ' opt_status ' 
                     + str(themodel.status) + ' runtime ' 
                     + str(all_data['runtime']) + ' iterations ' 
                     + str(all_data['round']) + '\n')

    summary_ws.close()

    log.joint(' solver runtime current round = %g\n' % (t1_solve - t0_solve) )
    log.joint(' overall time = %g\n' % all_data['runtime'] )
    log.joint(' bye.\n')
    return None

    
  elif ( themodel.status != GRB.status.OPTIMAL and 
         themodel.status != GRB.status.SUBOPTIMAL and
         themodel.status != GRB.status.NUMERIC):

    log.joint(' -> solver terminated with status ' + str(themodel.status) + '\n')

    all_data['runtime'] = time.time() - all_data['T0']

    log.joint(' writing casename, opt status and runtime to summary_ws.log\n')

    summary_ws = open("summary_ws.log","a+")
    summary_ws.write(' case ' + all_data['casename'] + ' opt_status ' 
                     + str(themodel.status) + ' runtime ' 
                     + str(all_data['runtime']) + ' iterations ' 
                     + str(all_data['round']) + '\n')

    summary_ws.close()

    log.joint(' solver runtime current round = %g\n' % (t1_solve - t0_solve) )
    log.joint(' overall time = %g\n' % all_data['runtime'] )
    log.joint(' bye.\n')
    return None

  all_data['objval']                  = themodel.ObjVal
  all_data['optstatus']               = themodel.status
  all_data['solvertime']              = t1_solve - t0_solve
  all_data['cumulative_solver_time'] += (t1_solve - t0_solve)  


def getsol_ampl(log,all_data):

  casename    = all_data['casename']
  filename    = 'ACsol_'+ casename +'.txt'

  try:
    thefile   = open(filename, "r")
    lines     = thefile.readlines()
    lenlines  = len(lines)
    thefile.close()
  except:
    log.stateandquit("Cannot open file " + datafilename)
    sys.exit("failure")

  branches     = all_data['branches']
  buses        = all_data['buses']
  IDtoCountmap = all_data['IDtoCountmap'] 

  sol_obj        = 0
  sol_vm         = {}
  sol_angle      = {}
  sol_cvalues    = {}
  sol_svalues    = {}
  sol_Pfvalues   = {}
  sol_Ptvalues   = {}
  sol_Qfvalues   = {}
  sol_Qtvalues   = {}
  sol_GenPvalues = {}

  linenum = 0
  log.joint(' reading file\n')
  while linenum < lenlines:
    thisline = lines[linenum].split()
    if thisline[0] == 'objvalue':
      sol_obj              = float(thisline[1])
    elif thisline[0] == 'bus':
      buscount             = int(thisline[1])
      bus                  = buses[buscount]
      sol_vm[bus]          = float(thisline[3])
      #sol_angle[bus]       = float(thisline[5]) 
      sol_cvalues[bus]     = sol_vm[bus]**2
    elif thisline[0] == 'branch':
      branchcount          = int(thisline[1])
      branch               = branches[branchcount]
      sol_Pfvalues[branch] = float(thisline[7])
      sol_Ptvalues[branch] = float(thisline[9])
      sol_Qfvalues[branch] = float(thisline[11])
      sol_Qtvalues[branch] = float(thisline[13])
      sol_cvalues[branch]  = float(thisline[15])
      sol_svalues[branch]  = float(thisline[17])
    elif thisline[0] == 'genid':
      genid      = int(thisline[1])
      gen_nodeID = int(thisline[3])
      sol_GenPvalues[genid] = float(thisline[5])

    linenum += 1

  all_data['sol_vm']       = sol_vm
  all_data['sol_angle']    = sol_angle
  all_data['sol_cvalues']  = sol_cvalues
  all_data['sol_svalues']  = sol_svalues
  all_data['sol_Pfvalues'] = sol_Pfvalues
  all_data['sol_Ptvalues'] = sol_Ptvalues
  all_data['sol_Qfvalues'] = sol_Qfvalues
  all_data['sol_Qtvalues'] = sol_Qtvalues
  all_data['sol_GenPvalues'] = sol_GenPvalues

  log.joint(' AMPL solution loaded\n')


def writesol(log,all_data):
    
  branches     = all_data['branches']
  buses        = all_data['buses']
  gens         = all_data['gens']
  IDtoCountmap = all_data['IDtoCountmap']
  cvalues      = all_data['cvalues']
  svalues      = all_data['svalues']
  GenPvalues   = all_data['GenPvalues']
  GenQvalues   = all_data['GenQvalues']
  Pfvalues     = all_data['Pfvalues']
  Ptvalues     = all_data['Ptvalues']
  Qfvalues     = all_data['Qfvalues']
  Qtvalues     = all_data['Qtvalues']
  i2fvalues    = all_data['i2fvalues']

  filename    = 'CPsol_' + all_data['casename'] + '.txt'
  thefile     = open(filename,'w+')

  log.joint('\n writing solution to ' + filename + '\n')

  machinename = platform.node()
  now = time.time()
  solver_version = 'Gurobi 11.0.3'  
  opsystem = "{} {} ({})".format(platform.system(), platform.release(), platform.version())
  processor = platform.processor()
  physical_cores = psutil.cpu_count(logical=False)
  logical_processors = psutil.cpu_count(logical=True)
  cores = f"{physical_cores} physical cores, {logical_processors} logical processors"
  ram = f"{round(psutil.virtual_memory().total / (1024 ** 3))} GB RAM"

  thefile.write('/CUTPLANEsolution : ' + all_data['casename'] + '\n')
  thefile.write('/Date : ' + str(time.strftime('%m-%d-%Y %H:%M:%S %Z', time.localtime(now))) + '\n')
  thefile.write('/MachineName : ' + machinename + '\n')
  thefile.write('/Processor : ' + processor + '\n')
  thefile.write('/OS : ' + opsystem + '\n')
  thefile.write('/Cores : ' + cores + '\n')
  thefile.write('/RAM : ' + ram + '\n')
  thefile.write('/Solver : ' + solver_version + '\n')
  thefile.write('objvalue ' + str(all_data['objval']) + '\n')
  thefile.write('round' + str(all_data['round']) + '\n')
    
  thefile.write('voltages:\n')

  for buscount in buses.keys():
    bus   = buses[buscount]
    v2val = cvalues[bus]
    vval  = (v2val)**(0.5)
    line  = 'bus ' + str(buscount) + ' M ' + str(vval) + '\n'
    thefile.write(line)

  thefile.write('power flows and cs variables:\n')

  for branchid in branches.keys():

    branch     = branches[branchid]
    f          = branch.f
    t          = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    Pfval      = Pfvalues[branch]
    Ptval      = Ptvalues[branch]
    Qfval      = Qfvalues[branch]
    Qtval      = Qtvalues[branch]
    cftval     = cvalues[branch]
    sftval     = svalues[branch]
    i2fval     = i2fvalues[branch]

    line = 'branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' Pft ' + str(Pfval) + ' Ptf ' + str(Ptval) + ' Qft ' + str(Qfval) + ' Qtf ' + str(Qtval) + ' cft ' + str(cftval) + ' sft ' + str(sftval) + ' i2ft ' + str(i2fval) + '\n'
    thefile.write(line)

  thefile.write('generation:\n')

  for genid in gens.keys():
    gen     = gens[genid] 
    nodeID  = gen.nodeID
    GenPval = GenPvalues[gen]
    GenQval = GenQvalues[gen]
    line_gen = 'genid ' + str(genid) + ' bus ' + str(nodeID) + ' GP ' + str(GenPval) + ' GQ ' + str(GenQval) + '\n'
    thefile.write(line_gen)

  thefile.close()

  log.joint(' done writing CUTPLANE solution to .txt file\n\n')



def writesol_allvars(log,all_data):

  branches     = all_data['branches']
  buses        = all_data['buses']
  gens         = all_data['gens']
  IDtoCountmap = all_data['IDtoCountmap']
  cvalues      = all_data['cvalues']
  svalues      = all_data['svalues']
  GenPvalues   = all_data['GenPvalues']
  GenQvalues   = all_data['GenQvalues']
  Pfvalues     = all_data['Pfvalues']
  Ptvalues     = all_data['Ptvalues']
  Qfvalues     = all_data['Qfvalues']
  Qtvalues     = all_data['Qtvalues']
  i2fvalues    = all_data['i2fvalues']


  filename    = 'CPsol_' + all_data['casename'] + '.sol'
  thefilevars = open(filename,'w+')

  log.joint('\n writing solution to ' + filename + '\n')

  machinename = platform.node()
  now = time.time()
  solver_version = 'Gurobi 11.0.3' 
  opsystem = "{} {} ({})".format(platform.system(), platform.release(), platform.version())
  processor = platform.processor()
  physical_cores = psutil.cpu_count(logical=False)
  logical_processors = psutil.cpu_count(logical=True)
  cores = f"{physical_cores} physical cores, {logical_processors} logical processors"
  ram = f"{round(psutil.virtual_memory().total / (1024 ** 3))} GB RAM"

  thefilevars.write('/CUTPLANEsolution : ' + all_data['casename'] + '\n')
  thefilevars.write('/Date : ' + str(time.strftime('%m-%d-%Y %H:%M:%S %Z', time.localtime(now))) + '\n')
  thefilevars.write('/MachineName : ' + machinename + '\n')
  thefilevars.write('/Processor : ' + processor + '\n')
  thefilevars.write('/OS : ' + opsystem + '\n')
  thefilevars.write('/Cores : ' + cores + '\n')
  thefilevars.write('/RAM : ' + ram + '\n')
  thefilevars.write('/Solver : ' + solver_version + '\n')
  thefilevars.write('/Objvalue ' + str(all_data['objval']) + '\n')

  for buscount in buses.keys():
    bus        = buses[buscount]
    f          = bus.nodeID
    v2value    = cvalues[bus]        

    v2name     = 'c_' + str(f) + '_' + str(f)
    v2line     = v2name + ' = ' + str(v2value) + '\n'
    thefilevars.write(v2line)

    IPvalue = - bus.Pd
    IQvalue = - bus.Qd
    IPname  = 'IP_' + str(bus.nodeID)
    IQname  = 'IQ_' + str(bus.nodeID)

    for gencounter in bus.genidsbycount:
      if gens[gencounter].status:
        gen      = gens[gencounter]
        IPvalue += GenPvalues[gen]
        IQvalue += GenQvalues[gen]

    IPline = IPname + ' = ' + str(IPvalue) + '\n'
    thefilevars.write(IPline)

    IQline = IQname + ' = ' + str(IQvalue) + '\n'
    thefilevars.write(IQline)


  for branchid in branches.keys():

    branch     = branches[branchid]
    f          = branch.f
    t          = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    Pfval      = Pfvalues[branch]
    Ptval      = Ptvalues[branch]
    Qfval      = Qfvalues[branch]
    Qtval      = Qtvalues[branch]
    cftval     = cvalues[branch]
    sftval     = svalues[branch]
    i2fval     = i2fvalues[branch]

    Pfname  = 'P_' + str(branchid) + '_' + str(f) + '_' + str(t)
    Ptname  = 'P_' + str(branchid) + '_' + str(t) + '_' + str(f)
    Qfname  = 'Q_' + str(branchid) + '_' + str(f) + '_' + str(t)
    Qtname  = 'Q_' + str(branchid) + '_' + str(t) + '_' + str(f)
    cftname = 'c_' + str(branchid) + '_' + str(f) + '_' + str(t)
    sftname = 's_' + str(branchid) + '_' + str(f) + '_' + str(t)
    i2fname = 'i2_'+ str(branchid) + '_' + str(f) + '_' + str(t)            

    Pfline  = Pfname + ' = ' + str(Pfval) + '\n'
    thefilevars.write(Pfline)

    Ptline  = Ptname + ' = ' + str(Ptval) + '\n'
    thefilevars.write(Ptline)

    Qfline  = Qfname + ' = ' + str(Qfval) + '\n'
    thefilevars.write(Qfline)

    Qtline  = Qtname + ' = ' + str(Qtval) + '\n'
    thefilevars.write(Qtline)

    cftline  = cftname + ' = ' + str(cftval) + '\n'
    thefilevars.write(cftline)

    sftline  = sftname + ' = ' + str(sftval) + '\n'
    thefilevars.write(sftline)

    i2fline  = i2fname + ' = ' + str(i2fval) + '\n'
    thefilevars.write(i2fline)                        


  for genid in gens.keys():
    gen     = gens[genid] 
    nodeID  = gen.nodeID
    GenPval = GenPvalues[gen]
    GenQval = GenQvalues[gen]
    GPname  = "GP_" + str(genid) + "_" + str(nodeID)
    GQname  = "GQ_" + str(genid) + "_" + str(nodeID)

    GPline  = GPname + ' = ' + str(GenPval) + '\n'
    thefilevars.write(GPline)

    GQline  = GQname + ' = ' + str(GenQval) + '\n'
    thefilevars.write(GQline)


  log.joint(' done writing CUTPLANE allvars solution to .sol file\n\n')

  thefilevars.close()
