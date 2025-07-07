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
from log import danoLogger
from gurobipy import *
import numpy as np
from myutils import *
import reader
import time
import math
from cuts_mtp import *
import os
import platform
import psutil

# This is the main function which starts the cutting-plane procedure

def gocutplane(log, all_data):

  #################### LOADING CASE PARAMETERS ################################
  
  formulation_start = time.time()
  themodel          = Model("Cutplane")
  buses             = all_data['buses']
  numbuses          = all_data['numbuses']
  branches          = all_data['branches']
  numbranches       = all_data['numbranches']
  gens              = all_data['gens']
  IDtoCountmap      = all_data['IDtoCountmap']
  T                 = all_data['T']
  casename          = all_data['casename']
  casetype          = all_data['casetype']
  loadsfilename     = all_data['loadsfilename']
  rampfilename      = all_data['rampfilename']
  
  ############################ LOAD SOLUTION ##################################

  if all_data['ampl_sol']:
    getsol_ampl_mtp(log,all_data)

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

  # Loading file with multi-period loads
  log.joint('Name of mtp file ' + loadsfilename + '\n')
  try:
    loads = open(loadsfilename,"r")
    Pd    = getloads(log,all_data,loads)
    all_data['Pd'] = Pd
    log.joint(" Loads obtained\n")
  except:
    log.joint(" File with mtp loads could not be found in '../data/mtploads/'")
    log.joint(" Please provide it\n")
    exit(0)
  
  # Loading file with ramping rates loads
  log.joint('Name of ramprates file ' + rampfilename + '\n')
  try:
    rampr          = open(rampfilename,"r")
    rampru, ramprd = getrampr(log,all_data,rampr)
    all_data['rampru'] = rampru
    all_data['ramprd'] = ramprd  
    log.joint(" Ramp rates obtained\n")
  except:
    log.joint(" File with ramping rates could not be found in '../data/ramprates/'")
    log.joint(" Please provide it\n")
    exit(0)
 
  
  for k in range(T):
    cvar[k]    = {}
    svar[k]    = {}
    Pvar_f[k]  = {}
    Qvar_f[k]  = {}
    Pvar_t[k]  = {}
    Qvar_t[k]  = {}
    Pinjvar[k] = {}
    Qinjvar[k] = {}
    GenPvar[k] = {}
    GenQvar[k] = {}
    GenTvar[k] = {}


  log.joint(' %d Time periods\n' %T)      
  log.joint(' Creating variables...\n')

  varcount = 0

  # Bus-related variables 
  for bus in buses.values():

    maxprod = bus.Vmax*bus.Vmax
    minprod = bus.Vmin*bus.Vmin
          
    ubound = maxprod
    lbound = minprod

    for k in range(T):
      Pubound, Plbound, Qubound, Qlbound = computebalbounds(log,all_data,bus,k)
      
      # cvar[k][bus]: represents the square of the voltage magnitude at bus 'bus'
      # in period 'k'
      # Pinjvar[k][bus]: represents the active power injection at bus 'bus' in 
      # period 'k'
      # Qinjvar[k][bus]: represents the reactive power injection at bus 'bus' in 
      # period 'k'
      cvar[k][bus] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                     name = "c_" + str(bus.nodeID) + "_" 
                                     + str(bus.nodeID) + "_" + str(k))
      Pinjvar[k][bus] = themodel.addVar(obj = 0.0, lb = Plbound, ub = Pubound, 
                                        name = "IP_"+str(bus.nodeID) + "_" + str(k))
      Qinjvar[k][bus] = themodel.addVar(obj = 0.0, lb = Qlbound, ub = Qubound, 
                                        name = "IQ_"+str(bus.nodeID) + "_" + str(k))
    
      varcount += 3

    # This loop goes over all generators (by genid) at bus 'bus'. 'genid' is a unique 
    # identifier for each generator
    for genid in bus.genidsbycount:
      gen = gens[genid]
      lower = gen.Pmin*gen.status
      upper = gen.Pmax*gen.status

      for k in range(T):

        # GenPvar[k][gen]: represents active power generated at generator
        # 'gen' in period 'k'
        GenPvar[k][gen] = themodel.addVar(obj = 0.0, lb = lower, ub = upper, 
                                          name = "GP_" + str(gen.count) + "_" 
                                          + str(gen.nodeID) + "_" + str(k))
        varcount += 1
      
      lower = gen.Qmin*gen.status
      upper = gen.Qmax*gen.status

      if bus.nodetype == 3:
        upper = GRB.INFINITY
        lower = -GRB.INFINITY

      for k in range(T):
        # GenQvar[k][gen]: represents reactive power generated at generator
        # 'gen' in period 'k'
        GenQvar[k][gen] = themodel.addVar(obj = 0.0, lb = lower, ub = upper, 
                                          name = "GQ_" + str(gen.count) + "_" 
                                          + str(gen.nodeID) + "_" + str(k))
        varcount += 1

  # Branch-related variables
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

    # Bounds for c variables                                                                                                                                          
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

    for k in range(T):
      # cvar[k][branch]: represents the linearization of the product 
      # v_i * v_j * cos(theta_i - theta_j), where 'i' denotes the 'from' bus
      # and 'j' denotes the 'to' bus, for branch 'branch' in period 'k', 
      # c.f. equation (3) in [1].

      cvar[k][branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                        name = "c_" + str(branchcount) + "_" 
                                        + str(f) + "_" + str(t) + "_" + str(k))
      varcount += 1
      
    # s variables
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

    for k in range(T):
      # svar[k][branch]: represents the linearization of the product 
      # v_i * v_j * sin(theta_i - theta_j), where 'i' denotes the 'from' bus
      # and 'j' denotes the 'to' bus, for branch 'branch' in period 'k', 
      # c.f. equation (3) in [1].

      svar[k][branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                        name = "s_" + str(branchcount) + "_" 
                                        + str(f) + "_" + str(t) + "_" + str(k))

      varcount += 1
    
  # Flow variables
  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]

    ubound = branch.limit 
    lbound = -branch.limit

    for k in range(T):
      # Pvar_f[k][branch]: given a branch 'branch' with 'from' bus 'f' and 'to'
      # bus 't', this variable represents the active power flow with 
      # orientation (f,t) in period 'k'
      # Qvar_f[k][branch]: given a branch 'branch' with 'from' bus 'f' and 'to'
      # bus 't', this variable represents the reactive power flow with 
      # orientation (f,t) in period 'k'
      # Similarly, variables Pvar_t[k][branch] and Qvar_t[k][branch] represent
      # active and reactive power flows, resp., with orientation (t,f) 
      # in period 'k'
      Pvar_f[k][branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                          name = "P_" + str(branch.count) + "_" 
                                          + str(f) + "_" + str(t) + "_" + str(k))
      Pvar_t[k][branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                          name = "P_" + str(branch.count) + "_" 
                                          + str(t) + "_" + str(f) + "_" + str(k))
      Qvar_f[k][branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                          name = "Q_" + str(branch.count) + "_" 
                                          + str(f) + "_" + str(t) + "_" + str(k))
      Qvar_t[k][branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                          name = "Q_" + str(branch.count) + "_" 
                                          + str(t) + "_" + str(f) + "_" + str(k))

      varcount +=4 

      
  # i2 variables
  if all_data['i2']:
    i2var_f = {}

    # This dictionary will contain the coeffs of |v|^2 associated to the
    # 'from' bus in the i2 definition, c.f. equation (10) in [1]. This
    # dictionary will be used for our heuristic described in section 4.5 in [1]
    alphadic = all_data['alphadic'] = {}
    
    for k in range(T):
      i2var_f[k] = {}
      

    for branch in branches.values():
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

      alpha            = ( g*g + b*b + bshunt * (b + (bshunt/4)) ) / (ratio**4)
      alphadic[branch] = alpha

      # By our heuristic in section 4.5 in [1], we explicitly define variable 
      # i2 for branch 'branch' if its 'alpha' coeff in equation (18) in [1]
      # is less or equal than 'rho_threshold'

      if alpha < all_data['rho_threshold']:
        bus_f        = buses[count_of_f]
        upperbound_f = branch.limit**2 / (bus_f.Vmin * bus_f.Vmin)
        
        for k in range(T):
          # i2var_f[k][branch]: represents the current squared at branch
          # 'branch' with orientation (f,t) at period 'k'
          i2var_f[k][branch] = themodel.addVar(obj = 0.0, lb = 0, ub = upperbound_f ,
                                               name = "i2_" + str(branch.count) + "_"
                                               + str(f) + "_" + str(t) + "_" + str(k))    
          varcount += 1
      
  themodel.update()
  log.joint('   %d variables added\n' %varcount)

  all_data['themodel']    = themodel
  all_data['cvar']        = cvar
  all_data['svar']        = svar
  all_data['GenPvar']     = GenPvar
  all_data['GenTvar']     = GenTvar
  all_data['Pvar_f']      = Pvar_f
  all_data['Pvar_t']      = Pvar_t
  all_data['Qvar_f']      = Qvar_f
  all_data['Qvar_t']      = Qvar_t
  all_data['Pd']          = Pd

  if all_data['i2']:
    all_data['i2var_f']   = i2var_f

  ############################## OBJECTIVE ####################################

  # Refer to section 2.1 in [1]
  log.joint(' Creating the objective...\n')
  varobjcount = 0

  # Constant term (fixed costs)
  constexpr   = LinExpr()
  constobjval = 0
  
  for gen in gens.values():  
    if gen.status > 0:
      constobjval += gen.costvector[gen.costdegree]
  
  constvar   = themodel.addVar(lb = 1.0, ub = 1.0,name = "constant")
  constexpr += T * constobjval * constvar
  
  varobjcount +=1

  # Linear and Quadratic terms
  lincostexpr = LinExpr()
  qcostexpr   = QuadExpr()
  
  for gen in gens.values():
    lincoeff = gen.costvector[gen.costdegree-1]
    
    for k in range(T):
      lincostexpr += lincoeff * GenPvar[k][gen]
      varobjcount +=1
      
      if gen.costdegree == 2 and gen.costvector[0] != 0:
        qcostexpr += gen.costvector[0]*GenPvar[k][gen]*GenPvar[k][gen]
        varobjcount +=1

  log.joint('   %d terms in the objective\n' %varobjcount)
  
  themodel.setObjective(constexpr + lincostexpr + qcostexpr)
  
  themodel.update()
  
  ############################# CONSTRAINTS ###################################

  log.joint(' Creating the constraints...\n')

  constrcount = 0
  count       = 0

  # Definition of Flow variables
  log.joint('  Active power flow variables definition\n')

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]

    for k in range(T):
      # Refer to equation (6a) in [1], 
      # Pft_k = Gff cff_k + Gft cft_k + Bft sft_k where 'k' denotes
      # the time-period
      constrname = "Pdef_"+str(branch.count)+"_"+str(f)+"_"+str(t)+"_"+str(k)
      expr = LinExpr()
      expr += branch.Gff*cvar[k][buses[count_of_f]]
      expr += branch.Gft*cvar[k][branch]
      expr += branch.Bft*svar[k][branch]
    
      themodel.addConstr(expr == Pvar_f[k][branch], name = constrname)

      # Refer to equation (6be) in [1], 
      # Ptf_k = Gtt ctt_k + Gtf cft_k + Btf stf_k 
      #       = Gtt ctt_k + Gtf cft_k - Btf sft_k
      # where 'k' denotes the time-period
      constrname = "Pdef_"+str(branch.count)+"_"+str(t)+"_"+str(f)+"_"+str(k)
      expr = LinExpr()
      expr += branch.Gtt*cvar[k][buses[count_of_t]]
      expr += branch.Gtf*cvar[k][branch]
      expr += -branch.Btf*svar[k][branch]

      themodel.addConstr(expr == Pvar_t[k][branch], name = constrname)
    
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

    for k in range(T):
      # Refer to equation (6c) in [1], 
      # Qft_k = -Bff cff_k - Bft cft_k + Gft sft_k where 'k' denotes
      # the time-period
      constrname = "Qdef_"+str(branch.count)+"_"+str(f)+"_"+str(t)+"_"+str(k)
      expr = LinExpr()
      expr += -branch.Bff*cvar[k][buses[count_of_f]]
      expr += -branch.Bft*cvar[k][branch]
      expr += +branch.Gft*svar[k][branch]

      themodel.addConstr(expr == Qvar_f[k][branch], name = constrname)

      # Refer to equation (6d) in [1], 
      # Qtf_k = -Btt ctt_k - Btf cft_k + Gtf stf_k
      #       = -Btt ctt_k - Btf cft_k - Gtf sft_k 
      # where 'k' denotes the time-period
      constrname = "Qdef_"+str(branch.count)+"_"+str(t)+"_"+str(f)+"_"+str(k)
      expr = LinExpr()
      expr += -branch.Btt*cvar[k][buses[count_of_t]]
      expr += -branch.Btf*cvar[k][branch]
      expr += -branch.Gtf*svar[k][branch]

      themodel.addConstr(expr == Qvar_t[k][branch], name = constrname)

      constrcount += 2
      count       += 2
  log.joint('   %d reactive power flow definition constraints added\n'%count)

  # Flow balance constraints
  log.joint('  active power injection constraints\n')
  count = 0

  for bus in buses.values():
    for k in range(T):
      # Refer to equation (1b) in [1]
      constrname = "PBaldef"+str(bus.nodeID)+"_"+str(k)
      expr = LinExpr()
      for branchid in bus.frombranchids.values():
        expr += Pvar_f[k][branches[branchid]]
      for branchid in bus.tobranchids.values():
        expr += Pvar_t[k][branches[branchid]]

      if ( (bus.Gs != 0) and ( ( len(bus.frombranchids) != 0 ) 
                               or ( len(bus.tobranchids) != 0 ) ) ):
        expr += bus.Gs*cvar[k][bus]

      themodel.addConstr(expr == Pinjvar[k][bus], name = constrname)

      constrcount += 1
      count       += 1
  
  log.joint('   %d active power injection constraints added\n'%count)
  log.joint('  reactive power injection constraints\n')
  count = 0

  for bus in buses.values():
    for k in range(T):
      # Refer to equation (1c) in [1]
      constrname = "QBaldef"+str(bus.nodeID)+"_"+str(k)
      expr = LinExpr()
      for branchid in bus.frombranchids.values():
        expr += Qvar_f[k][branches[branchid]]
      for branchid in bus.tobranchids.values():
        expr += Qvar_t[k][branches[branchid]]
 
      if ( (bus.Bs != 0) and ( ( len(bus.frombranchids) != 0 ) 
                               or ( len(bus.tobranchids) != 0 ) ) ):
        expr += (-bus.Bs)*cvar[k][bus]

      themodel.addConstr(expr == Qinjvar[k][bus], name = constrname)
    
      constrcount += 1
      count       += 1
      
  log.joint('   %d reactive power injection constraints added\n'%count)

  # Definition of Bus-injection variables
  log.joint('  Adding injection definition constraints...\n')
  count = 0

  for bus in buses.values():
    for k in range(T):
      # This constraint defines variable 'Pinjvar[k][bus]'
      # which represents net active power injection, i.e.,
      # total active power generation at bus 'bus' minus 
      # total load at bus 'bus', in period 'k'. 
      # See RHS of equation (1b) in [1]
      constrname = "Bus_PInj_"+str(bus.nodeID)+"_"+str(k)
      expr = LinExpr()

      if len(bus.genidsbycount) > 0:
        for genid in bus.genidsbycount:
          gen = gens[genid]
          expr += GenPvar[k][gen]

      themodel.addConstr(Pinjvar[k][bus] == expr - Pd[k][bus], name = constrname)

      # This constraint defines variable 'Qinjvar[k][bus]'
      # which represents net reactive power injection, i.e.,
      # total reactive power generation at bus 'bus' minus 
      # total load at bus 'bus', in period 'k'. 
      # See RHS of equation (1c) in [1]
      constrname = "Bus_QInj_"+str(bus.nodeID)+"_"+str(k)
      expr = LinExpr()

      if len(bus.genidsbycount) > 0:
        for genid in bus.genidsbycount:
          gen = gens[genid]
          expr += GenQvar[k][gen]

      themodel.addConstr(Qinjvar[k][bus] == expr - bus.Qd, name = constrname)

      constrcount += 2
      count       += 2

  log.joint('   %d power injection definitions added\n'%count)

  # Ramping up and down constraints
  log.joint('  Adding ramping constraints...\n')

  # Refer to equations (2a) and (2b) in [1]
  count   = 0
  abs_gen = {}
  for k in range(T-1):
    abs_gen[k] = {}

  for gen in gens.values():

    genid  = gen.count
    nodeID = gen.nodeID
    for k in range(T-1):
      constrname_rup   = "rup_" + str(genid) + "_" + str(nodeID) + "_" + str(k) + "_" + str(k+1)
      constrname_rdown = "rdown_" + str(genid) + "_" + str(nodeID) + "_" + str(k) + "_" + str(k+1)      
      rpu        = rampru[k][gen]
      rpd        = ramprd[k][gen]
      
      abs_gen[k][gen] = themodel.addVar(lb = 0, name = 'abs_gen_'+str(genid)+'_'+str(k))
      themodel.addConstr(GenPvar[k+1][gen] - GenPvar[k][gen] - rpu * abs_gen[k][gen] <= 0, name = constrname_rup)
      themodel.addConstr(GenPvar[k][gen] - rpd * abs_gen[k][gen] - GenPvar[k+1][gen] <= 0, name = constrname_rdown)

      themodel.addConstr(GenPvar[k][gen] - abs_gen[k][gen] <= 0, name = constrname_rup + '_1')
      themodel.addConstr(- GenPvar[k][gen] - abs_gen[k][gen] <= 0, name = constrname_rup + '_2')      

      count += 4
  
  # Definition i2 variables
  if all_data['i2']:
    constrcount += i2_def(log,all_data)

  # Jabr inequalities
  # Refer to equation (1k) in [1]
  if all_data['jabr_inequalities']:
    constrcount += jabr_inequalities(log,all_data)

  # i2 inequalities
  # Refer to equation (1k) in [1]
  if all_data['i2_inequalities']:
    constrcount += i2_inequalities(log,all_data)

  # Limit constraints
  # Refer to equation (1k) in [1]
  if all_data['limit_inequalities']:
    constrcount += limit_inequalities(log,all_data)
  
  log.joint('  %d constraints added\n'%constrcount)
    
  themodel.update()

  formulation_end = time.time()

  all_data['formulation_time'] = formulation_end - formulation_start
  all_data['numvars']          = themodel.NumVars
  all_data['numconstrs']       = themodel.NumConstrs
  
  log.joint(' Formulation time: %g\n' % all_data['formulation_time'])
  log.joint(' numvars ' + str(all_data['numvars']) + ' numconstrs ' + str(all_data['numconstrs']) + '\n')
  
  # Write model to a .lp file
  if all_data['writelps']:
    log.joint(' writing to lpfile ' + all_data['lpfilename'] + '\n')  
    themodel.write(all_data['lpfilename'])

  
  ###################### INIT DATA STRUCTURES FOR CUTS ########################

  # These dictionaries will collect information of all of the cuts
  # computed throughout our cutting-plane procedure

  if all_data['jabrcuts']:
    jabr_cuts_info = all_data['jabr_cuts_info']
    for k in range(T):
      jabr_cuts_info[k] = {}
      for branch in branches.values():
        jabr_cuts_info[k][branch] = {}

  if all_data['i2cuts']:
    i2_cuts_info = all_data['i2_cuts_info']
    for k in range(T):
      i2_cuts_info[k] = {}
      for branch in branches.values():
        i2_cuts_info[k][branch] = {}

  if all_data['limitcuts']:
    limit_cuts_info = all_data['limit_cuts_info']
    for k in range(T):
      limit_cuts_info[k] = {}
      for branch in branches.values():
        limit_cuts_info[k][branch] = {}

  ######################## FIXING/WRITING AN AC SOLUTION ######################

  # The following functions use ac AC solution previously loaded via 'ampl_sol'
  # fixflows: This function fixes the flows (active and reactive power) up to 
  # some given tolerance using an AC solution 
  # fixcs: fixes c and s variables using an AC solution 

  if all_data['fixflows']:
    fixflows(log,all_data)
    if all_data['fixcs'] == 0:
      return None

  if all_data['fixcs']:
    fixcs(log,all_data)
    return None

  ########################### SOLVER PARAMETERS ###############################

  # By default we run Gurobi with the barrier algorithm and crossover disabled

  themodel.Params.Method    = all_data['solver_method']
  themodel.Params.Crossover = all_data['crossover'] 
  themodel.Params.LogFile   = all_data['mylogfile']
  themodel.Params.TimeLimit = all_data['max_time']

  if all_data['solver_method'] == 2:
    themodel.Params.BarHomogeneous = 1
    themodel.Params.BarConvTol     = all_data['barconvtol']
    themodel.Params.FeasibilityTol = all_data['feastol']
    themodel.Params.OptimalityTol  = all_data['opttol']
    
  themodel.Params.NumericFocus = 1
  themodel.Params.OutPutFlag = 1
  
  ######################### READING AND LOADING CUTS ##########################

  # This procedure adds previously computed cuts to the current optimization
  # instance. The function 'add_cuts_ws' is used if multiple cutting-plane
  # rounds want to be run after loading the cuts, and 'add_cuts' if only one 
  # iteration is needed

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
      log.joint(' model with precomputed written to .lp file\n\n')
  
  ########################## CUTPLANE MAIN LOOP ###############################

  all_data['round']                  = 1
  all_data['runtime']                = time.time() - all_data['T0']
  all_data['round_time']             = time.time()
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

    log.joint(' Storing current solution ...\n')

    all_data['Pfvalues']   = {}
    all_data['Qfvalues']   = {}
    all_data['Ptvalues']   = {}
    all_data['Qtvalues']   = {}
    all_data['cvalues']    = {}
    all_data['svalues']    = {}
    all_data['GenPvalues'] = {}
    all_data['GenQvalues'] = {}
    
    for k in range(T):
      Pfvalues   = {}
      Qfvalues   = {}
      Ptvalues   = {}
      Qtvalues   = {}
      cvalues    = {}
      svalues    = {}
      GenPvalues = {}
      GenQvalues = {}

      for bus in buses.values():
        cvalues[bus] = cvar[k][bus].X
      
      for branch in branches.values():
        Pfvalues[branch]   = Pvar_f[k][branch].X
        Qfvalues[branch]   = Qvar_f[k][branch].X
        Ptvalues[branch]   = Pvar_t[k][branch].X
        Qtvalues[branch]   = Qvar_t[k][branch].X
        cvalues[branch]    = cvar[k][branch].X
        svalues[branch]    = svar[k][branch].X

      for gen in gens.values():
        GenPvalues[gen] = GenPvar[k][gen].X
        GenQvalues[gen] = GenQvar[k][gen].X        

      all_data['Pfvalues'][k]   = Pfvalues
      all_data['Qfvalues'][k]   = Qfvalues
      all_data['Ptvalues'][k]   = Ptvalues
      all_data['Qtvalues'][k]   = Qtvalues
      all_data['cvalues'][k]    = cvalues
      all_data['svalues'][k]    = svalues
      all_data['GenPvalues'][k] = GenPvalues
      all_data['GenQvalues'][k] = GenQvalues      
      
    
    if all_data['i2']:
      all_data['i2fvalues'] = {}
      for k in range(T):
        i2fvalues = {}
        for branch in branches.values():   # check, loop over selected br
          alpha = all_data['alphadic'][branch]
          if alpha < all_data['rho_threshold']:
            i2fvalues[branch] = i2var_f[k][branch].X
        all_data['i2fvalues'][k] = i2fvalues
        
    log.joint(' done storing values\n')
     
    ########################## CHECK OBJ IMPROVEMENT ##########################

    # Here we check the objective had a relative improvement (wrt the 
    # previous objective) of at least 'ftol'. If this is not the case
    # then we increase the 'ftol_counter', else we reset it to 0

    if ((all_data['objval'] - oldobj)/oldobj) < all_data['ftol']:
      all_data['ftol_counter'] += 1
    else:
      all_data['ftol_counter'] = 0

    oldobj              = all_data['objval']
    all_data['runtime'] = time.time() - all_data['T0']

    ########################### ROUND STATISTICS ##############################

    cutplane_stats(log,all_data)

    ######################### SUMMARY EXPERIMENTS #############################

    # Here we log some statistics of the current round

    all_data['runtime'] = time.time() - all_data['T0']

    log.joint("\n writing casename, opt stauts, obj and " +
              "runtime to summary_ws.log\n")

    summary_ws   = open("summary_ws.log","a+")

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

    ############################ GET DUALS #################################

    # This function gets the dual variables associated to the active power 
    # balance constraints

    if all_data['getduals'] and (themodel.status != GRB.status.NUMERIC):
      getduals(log,all_data)
    
    ############################ TERMINATION #################################

    if (all_data['round'] >= all_data['max_rounds']):

      writesol_and_lps(log,all_data)

      summary_ws = open("summary_ws.log","a+") 
      summary_ws.write(' rounds limit reached!\n\n')
      summary_ws.close()
      log.joint(' rounds limit reached!\n')
      log.joint(' bye\n')
      return 0
          
    if all_data['runtime'] > all_data['max_time']:

      writesol_and_lps(log,all_data)

      summary_ws = open("summary_ws.log","a+")
      summary_ws.write(' time limit reached!\n\n')
      summary_ws.close()
      log.joint(' time limit reached!\n')
      log.joint(' bye\n')
      return 0

    if (all_data['ftol_counter'] > all_data['ftol_iterates']):

      writesol_and_lps(log,all_data)
     
      summary_ws = open("summary_ws.log","a+") 
      summary_ws.write(' poor consecutive obj improvement limit reached!\n\n')
      summary_ws.close()
      log.joint(' poor consecutive obj improvement limit reached\n')
      log.joint(' bye\n')
      return 0

    ############################### CUTS ######################################

    # Cut computations and management
    cutplane_cuts(log,all_data)

    # Cut statistics
    cutplane_cutstats(log,all_data)
    
    themodel.update()

    log.joint(' model updated\n')
    log.joint('\n')

    ############################### WRITE CUTS ################################

    # This function writes all the current cuts to a .txt file
    if all_data['writecuts']:
      write_cuts(log,all_data)

    ############################### WRITE LPS #################################
    
    # If this parameter is turned on then we write to a .lp file our 
    # current relaxation

    if all_data['writelps']:
      name = 'post_cuts' + '_' + str(all_data['round']) + '.lp'
      themodel.write(name)
      log.joint(' model with new cuts written to .lp file\n')

        
    ###########################################################################
                                              
    all_data['round']      += 1
    all_data['round_time']  = time.time()


###############################################################################

# Other Functions

# Writes down solutions, retrieves duals variables of active power balance 
# constraints, and writes to an .lp file our last linearly-constrained 
# relaxation

def writesol_and_lps(log,all_data):

  themodel = all_data['themodel']
  casename = all_data['casename']
  casetype = all_data['casetype']
  T        = all_data['T']

  # We write down our current solution to two files: the first
  # function creates a readable .txt where variables are sorted 
  # by type and index (i.e., voltages, power flows, generation); 
  # while the second function creates a .sol file where variables
  # are written (in arbitrary order) in the format: 
  # 'variable name' = 'variable value'
  if all_data['writesol']:
    writesol(log,all_data)
    writesol_allvars(log,all_data)

  # We print the duals associated to the active power balance
  # constraints and we write them down to a table
  if all_data['getduals'] and (themodel.status != GRB.status.NUMERIC):
    print_duals(log,all_data)
    print_duals_table(log,all_data)
  
  # If this parameter is turned on, then we write to an .lp file
  # our last linearly-constrained relaxation
  if all_data['writelastLP']:
    log.joint(' writing down last lp...\n')        
    themodel.write(casename + '_' + str(T) + '_' + casetype + "_last.lp")

  return 0


# Fixes the flows of some previously computed AC solution up to some
# given tolerance

def fixflows(log,all_data):

  if (all_data['ampl_sol'] == 0):
    log.joint(' cannot fix flows since no solution has been loaded!\n')
    return None

  themodel    = all_data['themodel']
  branches    = all_data['branches']
  T           = all_data['T']
  tolerance   = all_data['fix_tolerance']
  Pvar_f      = all_data['Pvar_f']
  Pvar_t      = all_data['Pvar_t']
  Qvar_f      = all_data['Qvar_f']
  Qvar_t      = all_data['Qvar_t']

  sol_Pfvalues = all_data['sol_Pfvalues']
  sol_Ptvalues = all_data['sol_Ptvalues']
  sol_Qfvalues = all_data['sol_Qfvalues']
  sol_Qtvalues = all_data['sol_Qtvalues']

  for branch in branches.values():
    for k in range(T):
      sol_Pf = sol_Pfvalues[k][branch]
      sol_Pt = sol_Ptvalues[k][branch]
      sol_Qf = sol_Qfvalues[k][branch]
      sol_Qt = sol_Qtvalues[k][branch]

      #Pf
      ubound_Pf = sol_Pf + tolerance
      lbound_Pf = sol_Pf - tolerance

      Pvar_f[k][branch].setAttr("ub",ubound_Pf)
      Pvar_f[k][branch].setAttr("lb",lbound_Pf)

      #Pt
      ubound_Pt = sol_Pt + tolerance
      lbound_Pt = sol_Pt - tolerance

      Pvar_t[k][branch].setAttr("ub",ubound_Pt)
      Pvar_t[k][branch].setAttr("lb",lbound_Pt)

      #Qf
      ubound_Qf = sol_Qf + tolerance
      lbound_Qf = sol_Qf - tolerance

      Qvar_f[k][branch].setAttr("ub",ubound_Qf)
      Qvar_f[k][branch].setAttr("lb",lbound_Qf)

      #Qt
      ubound_Qt = sol_Qt + tolerance
      lbound_Qt = sol_Qt - tolerance

      Qvar_t[k][branch].setAttr("ub",ubound_Qt)
      Qvar_t[k][branch].setAttr("lb",lbound_Qt)

  themodel.update()
  themodel.write('fixflows.lp')
  log.joint('check fixflows.lp\n')  

# Fixes the c and s variables of some previously computed AC solution
# up to some given tolerance

def fixcs(log,all_data):

  if (all_data['ampl_sol'] == 0):
    log.joint(' cannot fix flows since no solution has been loaded!\n')
    return None

  T           = all_data['T']
  themodel    = all_data['themodel']
  tolerance   = all_data['fix_tolerance']
  buses       = all_data['buses']
  branches    = all_data['branches']
  cvar        = all_data['cvar']
  svar        = all_data['svar']
  sol_cvalues = all_data['sol_cvalues']
  sol_svalues = all_data['sol_svalues']

  for bus in buses.values():

    for k in range(T):
      sol_v2 = all_data['sol_cvalues'][k][bus]
      ubound_v2 = sol_v2 + tolerance
      lbound_v2 = sol_v2 - tolerance

      cvar[k][bus].setAttr("ub",ubound_v2)
      cvar[k][bus].setAttr("lb",lbound_v2)

  for branch in branches.values():
    for k in range(T):
    
      sol_c = all_data['sol_cvalues'][k][branch]
      sol_s = all_data['sol_svalues'][k][branch]

      #c
      ubound_c = sol_c + tolerance
      lbound_c = sol_c - tolerance

      cvar[k][branch].setAttr("ub",ubound_c)
      cvar[k][branch].setAttr("lb",lbound_c)

      #s
      ubound_s = sol_s + tolerance
      lbound_s = sol_s - tolerance

      svar[k][branch].setAttr("ub",ubound_s)
      svar[k][branch].setAttr("lb",lbound_s)

  themodel.update()
  themodel.write('fixCS.lp')
  log.joint('check fixCS.lp\n')

# Computes bounds for active and reactive power injections 

def computebalbounds(log, all_data, bus, k):

  # We first get max/min generation
  baseMVA = all_data['baseMVA']
  gens = all_data['gens']
  Pd   = all_data['Pd']

  Pubound = Plbound = 0
  Qubound = Qlbound = 0

  for gencounter in bus.genidsbycount:
    if gens[gencounter].status:
      Pubound += gens[gencounter].Pmax
      Plbound += gens[gencounter].Pmin
      Qubound += gens[gencounter].Qmax
      Qlbound += gens[gencounter].Qmin

      if bus.nodetype == 3:
        Qubound = + GRB.INFINITY
        Qlbound = - GRB.INFINITY
        
  Pubound -= Pd[k][bus]
  Plbound -= Pd[k][bus]
  
  Qubound -= bus.Qd
  Qlbound -= bus.Qd

  if bus.nodetype == 4:
    Pubound = Plbound = Qubound = Qlbound = 0
    Pd[k][bus] = 0
      
  return Pubound, Plbound, Qubound, Qlbound


# Computes the value of the i2 variable of a given branch using squared 
# voltages of 'from' and 'to' buses (mp_cbusf,mp_cbust) and the corresponding 
# c and s values (mp_c,mp_s). See equation (29) in [1]

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


# Adds the Jabr rotated cone inequalities, see equation (5) in [1]

def jabr_inequalities(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  buses          = all_data['buses']
  cvar           = all_data['cvar']
  svar           = all_data['svar']
  IDtoCountmap   = all_data['IDtoCountmap']
  FeasibilityTol = all_data['tolerance']
  T              = all_data['T']
  
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
  maxk         = -1
  counter_jabr = 0

  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]

    if all_data['ampl_sol'] and all_data['jabr_validity']:
      for k in range(T):
        sol_c         = all_data['sol_cvalues'][k][branch]
        sol_s         = all_data['sol_svalues'][k][branch]
        sol_cbusf     = all_data['sol_cvalues'][k][buses[count_of_f]]
        sol_cbust     = all_data['sol_cvalues'][k][buses[count_of_t]]
        relviolation = violation = sol_c * sol_c + sol_s * sol_s - sol_cbusf * sol_cbust
        #relviolation = violation / ( ( coeff_Pft**2 + coeff_Qft**2 + coeff_cff**2 + coeff_i2ft**2 )**0.5 )
        if relviolation > maxviolation:
          maxviolation = relviolation
          maxbranch    = branch.count
          maxbusf      = f
          maxbust      = t
          maxk         = k

        if relviolation > FeasibilityTol:
          violated += 1
          log.joint('   WARNING, the Jabr-inequality associated to branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' is violated by the AC solution!\n')
          log.joint('   violation ' + str(violation) + '\n')
          log.joint('   relative violation ' + str(relviolation) + '\n')
          log.joint('   values (AC solution) ' + ' cft ' + str(sol_c) + ' sft ' + str(sol_s) + ' cff ' + str(sol_cbusf) + ' ctt ' + str(sol_cbust) + ' k ' + str(k) + '\n' )

        else:
          log.joint('   AC solution satisfies Jabr-inequality at branch ' + str(branchcount) + ' with slack ' + str(relviolation) + '\n')

    for k in range(T):
      trigexp       = QuadExpr()
      constrname    = "jabr_"+str(branchcount)+"_"+str(f)+"_"+str(t)+"_"+str(k)
      trigexp      += cvar[k][branch]*cvar[k][branch] + svar[k][branch]*svar[k][branch] - cvar[k][buses[count_of_f]]*cvar[k][buses[count_of_t]]

      themodel.addConstr(trigexp <= 0, name = constrname)
      counter_jabr += 1
    
  if all_data['ampl_sol'] and all_data['jabr_validity']:
    log.joint('  max violation of Jabr-inequalities by AC solution ' + str(maxviolation) + ' at branch ' + str(maxbranch) + ' f ' + str(maxbusf) + ' t ' + str(maxbust) + '\n')
    log.joint('  number of violated Jabr-inequalities ' + str(violated) + '\n')
    
  log.joint('   %d Jabr inequalities added\n'%counter_jabr) 

  return counter_jabr

# Defines the i2 variable, see equations (19a), (19b) and (29) in [1]

def i2_def(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  buses          = all_data['buses']
  i2var_f        = all_data['i2var_f']
  cvar           = all_data['cvar']
  svar           = all_data['svar']
  IDtoCountmap   = all_data['IDtoCountmap']
  T              = all_data['T']
  alphadic       = all_data['alphadic']
  
  counter_i2def  = 0
  counter_i2con  = 0
  
  log.joint('  i2 variables definition and i2 linear inequalities\n')

  for branch in branches.values():
    expr_f      = LinExpr()
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
    bus_f       = buses[count_of_f]
    bus_t       = buses[count_of_t]
    
    alpha = alphadic[branch]
    beta  = ( g*g + b*b ) / (ratio**2)
    gamma = ( math.cos(angle) * ( - 2 * (g*g + b*b) - b * bshunt ) + math.sin(angle) * ( - g * bshunt) ) / (ratio**3)
    zeta  = ( math.sin(angle) * ( - 2 * (g*g + b*b) - b * bshunt ) - math.cos(angle) * ( - g* bshunt) ) / (ratio**3)

    
    if alpha < all_data['rho_threshold']:
      for k in range(T):
        constrname_f = 'i2def_'+str(branch.count)+"_"+str(f) + "_" + str(t) + "_" + str(k)
        expr_f = alpha * cvar[k][bus_f] + beta * cvar[k][bus_t] + gamma * cvar[k][branch] + zeta * svar[k][branch]

        themodel.addConstr(expr_f == i2var_f[k][branch],name = constrname_f) 
        counter_i2def += 1
    else:
      upperbound_f = branch.limit**2 / (bus_f.Vmin * bus_f.Vmin)
      newalpha     = 1
      newbeta      = beta/alpha
      newgamma     = gamma/alpha
      newzeta      = zeta/alpha
      RHS          = upperbound_f / alpha
 
      for k in range(T):
        uppconstrname_f = 'uppi2_'+str(branchcount)+"_"+str(f) + "_" + str(t) + "_" + str(k)
        lowconstrname_f = 'lowi2_'+str(branchcount)+"_"+str(f) + "_" + str(t) + "_" + str(k)

        expr_f = cvar[k][bus_f] + newbeta * cvar[k][bus_t] + newgamma * cvar[k][branch] + newzeta * svar[k][branch]
        themodel.addConstr(expr_f <= RHS, name = uppconstrname_f)
        themodel.addConstr(expr_f >= 0,   name = lowconstrname_f)
        counter_i2con += 2
     
  log.joint('   %d i2 definition constraints added\n'%counter_i2def) 
  log.joint('   %d i2 linear constraints added\n'%counter_i2con)
  
  return counter_i2def + counter_i2con


# Adds the i2 rotated cone inequalities, see equation (9) in [1]

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
  T              = all_data['T']
  alphadic       = all_data['alphadic']

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
  maxk         = -1
  
  counter_i2   = 0

  for branch in branches.values():
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
    bus_f       = buses[count_of_f]
    bus_t       = buses[count_of_t]
    
    alpha = alphadic[branch]
    beta  = ( g*g + b*b ) / (ratio**2)
    gamma = ( math.cos(angle) * ( - 2 * (g*g + b*b) - b * bshunt ) + math.sin(angle) * ( - g * bshunt) ) / (ratio**3)
    zeta  = ( math.sin(angle) * ( - 2 * (g*g + b*b) - b * bshunt ) - math.cos(angle) * ( - g* bshunt) ) / (ratio**3)

    
    if alpha < all_data['rho_threshold']:
      if all_data['ampl_sol'] and all_data['i2_validity']:
        for k in range(T):
          sol_Pf        = all_data['sol_Pfvalues'][k][branch]
          sol_Qf        = all_data['sol_Qfvalues'][k][branch]
          sol_c         = all_data['sol_cvalues'][k][branch]
          sol_s         = all_data['sol_svalues'][k][branch]
          sol_cbusf     = all_data['sol_cvalues'][k][buses[count_of_f]]
          sol_cbust     = all_data['sol_cvalues'][k][buses[count_of_t]]
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
            maxk         = k

          if relviolation > FeasibilityTol:
            violated += 1
            log.joint('   WARNING, the i2 inequality associated to branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' is violated by the AC solution!\n')
            log.joint('   violation ' + str(violation) + '\n')
            log.joint('   relative violation ' + str(relviolation) + '\n')
            log.joint('   values (AC solution) ' + ' Pft ' + str(sol_Pf) + ' Qft ' + str(sol_Qf) + ' cff ' + str(sol_cbusf) + ' i2ft ' + str(sol_i2f) + ' k ' + str(k) + '\n' )
            #breakexit('check!')
          else:
            log.joint('   AC solution satisfies i2 inequality at branch ' + str(branchcount) + ' with slack ' + str(relviolation) + '\n')
            
      for k in range(T):
        counter_i2 += 1
        trigexp     = QuadExpr()
        constrname  = "i2_"+str(branchcount)+"_"+str(f)+"_"+str(t)+"_"+str(k)
        trigexp    += Pvar_f[k][branch]**2 + Qvar_f[k][branch]**2 - cvar[k][buses[count_of_f]] * i2var_f[k][branch]
        themodel.addConstr(trigexp <= 0, name = constrname)

  if all_data['ampl_sol'] and all_data['i2_validity']:
    log.joint('  max violation of i2 inequalities by AC solution ' + str(maxviolation) + ' at branch ' + str(maxbranch) + ' f ' + str(maxbusf) + ' t ' + str(maxbust) + ' k ' + str(maxk) + '\n')
    log.joint('  values (AC solution) ' + ' Pft ' + str(maxPf) + ' Qft ' + str(maxQf) + ' cff ' + str(maxcff) + ' i2ft ' + str(maxi2f) + '\n' )
    log.joint('  number of violated i2 inequalities ' + str(violated) + '\n')
    #breakexit('  check i2 violation')

  log.joint('   %d i2 inequalities added\n'%counter_i2) 

  return counter_i2  

# Adds the line limit inequalities, see equation (1k) in [1]

def limit_inequalities(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  Pvar_f         = all_data['Pvar_f']
  Pvar_t         = all_data['Pvar_t']
  Qvar_f         = all_data['Qvar_f']
  Qvar_t         = all_data['Qvar_t']
  T              = all_data['T']
  
  counter_limit  = 0

  log.joint('  limit inequalities\n')

  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t

    for k in range(T):
      constrname = "limit_f_"+str(branchcount)+"_"+str(f)+"_"+str(t)+"_"+str(k)
      limexp = QuadExpr()
      limexp += Pvar_f[k][branch]*Pvar_f[k][branch] + Qvar_f[k][branch]*Qvar_f[k][branch]
      themodel.addConstr(limexp <= branch.limit**2, name = constrname)

      constrname = "limit_t_"+str(branchcount)+"_"+str(t)+"_"+str(f)+"_"+str(k)
      limexp = QuadExpr()
      limexp += Pvar_t[k][branch]*Pvar_t[k][branch] + Qvar_t[k][branch]*Qvar_t[k][branch]
      themodel.addConstr(limexp <= branch.limit**2, name = constrname)

      counter_limit += 2

  log.joint('   %d limit inequalities added\n'%counter_limit) 

  return counter_limit

# Prints round statistics of the cutting-plane algorithm

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
  log.joint(' BarConvTol = ' + str(themodel.params.BarConvTol) 
            + ' FeasTol = ' + str(themodel.params.FeasibilityTol) 
            + ' OptTol = ' + str(themodel.params.OptimalityTol) + '\n') 

  if all_data['optstatus'] == 2 or all_data['optstatus'] == 13:
    log.joint(' max (unscaled/scaled) dual constraint error =\n'
              + '  ' + str(themodel.DualResidual) + ' / '
              + str(themodel.DualSResidual) + '\n')

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

    log.joint(' total number of cuts = ' + str(all_data['num_jabr_cuts']
                                               + all_data['num_i2_cuts']
                                               + all_data['num_limit_cuts'])
              + '\n')
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


# Prints cut statistics of the cutting-plane algorithm

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


  log.joint(' ---\n')
  log.joint(' total number of cuts = ' + str(all_data['num_jabr_cuts']
                                             + all_data['num_i2_cuts']
                                             + all_data['num_limit_cuts'])
            + '\n')
  log.joint(' ************************************************************\n\n') 


# Prints cut statistics of the current round of our cutting-plane procedure

def cutplane_cuts(log,all_data):

  log.joint('\n starting cut procedure ...\n')

  t0_cuts = time.time()

  t0_jabr = time.time()

  if all_data['jabrcuts']:
    jabr_cuts(log,all_data)

    if all_data['NO_jabrs_violated']:
      if all_data['threshold'] > all_data['tolerance']:
        all_data['threshold'] *= 1e-01
        log.joint(' threshold updated to ' + str(all_data['threshold'])
                  + '\n' )

  t1_jabr = time.time()

  log.joint(' time spent on Jabr-cuts ' + str(t1_jabr - t0_jabr) + '\n')

  t0_i2 = time.time()

  if all_data['i2cuts']:
    i2_cuts(log,all_data)

    if all_data['NO_i2_cuts_violated']:
      if all_data['threshold_i2'] > all_data['tolerance']:
        all_data['threshold_i2'] *= 1e-01
        log.joint(' threshold updated to ' + str(all_data['threshold_i2'])
                  + '\n' )

  t1_i2 = time.time()
  log.joint(' time spent on i2-cuts ' + str(t1_i2 - t0_i2) + '\n')


  t0_lim = time.time()

  if all_data['limitcuts']:
    limit_cuts(log,all_data)
    if all_data['NO_limit_cuts_violated']:
      if all_data['threshold_limit'] > all_data['tolerance']:
        all_data['threshold_limit'] *= 1e-01
        log.joint(' threshold updated to ' + str(all_data['threshold_limit'])
                  + '\n' )

  t1_lim = time.time()
  log.joint(' time spent on lim-cuts ' + str(t1_lim - t0_lim) + '\n')

  
  t1_cuts = time.time()

  log.joint('\n time spent on cuts ' + str(t1_cuts - t0_cuts) + '\n')


# Calls our cut management heuristics

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

  log.joint(' time spent on cut management ' + str(t1_cutmanag - t0_cutmanag)
            + '\n')


# Optimizes the current model

def cutplane_optimize(log,all_data):

  themodel = all_data['themodel']

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
    exit(0)

  elif themodel.status == GRB.status.INFEASIBLE:
    log.joint(' -> LP infeasible\n')

    all_data['runtime'] = time.time() - all_data['T0']

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
    exit(0)

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
    exit(0)

  all_data['objval']                  = themodel.ObjVal
  all_data['optstatus']               = themodel.status
  all_data['solvertime']              = t1_solve - t0_solve
  all_data['cumulative_solver_time'] += (t1_solve - t0_solve)
  if (themodel.status != GRB.status.NUMERIC):
    all_data['dinfs']                   = themodel.DualResidual
    all_data['dinfs_scaled']            = themodel.DualSResidual
  else:
    all_data['dinfs'] = -1
    all_data['dinfs_scaled'] = -1


# Loads a previously computed AC solution

def getsol_ampl_mtp(log,all_data):

  T           = all_data['T']
  casename    = all_data['casename']
  casetype    = all_data['casetype']

  filename      = 'ACsol_' + all_data['casename'] + '_' + str(T) + '_' + casetype + '.txt'
  
  try:
    thefile   = open(filename, "r")
    lines     = thefile.readlines()
    lenlines  = len(lines)
    thefile.close()
  except:
    log.stateandquit("cannot open file " + filename)
    sys.exit("failure")

  branches     = all_data['branches']
  buses        = all_data['buses']
  IDtoCountmap = all_data['IDtoCountmap'] 
  gens         = all_data['gens']
  tolerance    = 1e-05

  
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
  sol_GenQvalues = {}
  
  for k in range(T):
      sol_vm[k]         = {}
      sol_angle[k]      = {}
      sol_cvalues[k]    = {}
      sol_svalues[k]    = {}
      sol_Pfvalues[k]   = {}
      sol_Ptvalues[k]   = {}
      sol_Qfvalues[k]   = {}
      sol_Qtvalues[k]   = {}
      sol_GenPvalues[k] = {}
      sol_GenQvalues[k] = {}
      
  linenum = 0
  log.joint(' reading file\n')
  while linenum < lenlines:
    thisline = lines[linenum].split()
    if thisline[0] == 'objvalue':
      sol_obj              = float(thisline[1])
    elif thisline[0] == 'bus':
      buscount             = int(thisline[1])
      k                    = int(thisline[7])
      bus                  = buses[buscount]
      sol_vm[k][bus]       = float(thisline[3])
      sol_angle[k][bus]    = float(thisline[5]) 
      sol_cvalues[k][bus]  = float(thisline[3])**2
    elif thisline[0] == 'branch':   
      branchcount             = int(thisline[1])
      k                       = int(thisline[19])
      branch                  = branches[branchcount]
      sol_Pfvalues[k][branch] = float(thisline[7])
      sol_Ptvalues[k][branch] = float(thisline[9])
      sol_Qfvalues[k][branch] = float(thisline[11])
      sol_Qtvalues[k][branch] = float(thisline[13])
      sol_cvalues[k][branch]  = float(thisline[15])
      sol_svalues[k][branch]  = float(thisline[17])      
    elif thisline[0] == 'genid':
      genid      = int(thisline[1])
      k          = int(thisline[9])
      gen_nodeID = int(thisline[3])
      sol_GenPvalues[k][genid] = float(thisline[5])
      sol_GenQvalues[k][genid] = float(thisline[7])      

    linenum += 1

  all_data['sol_obj']        = sol_obj
  all_data['sol_vm']         = sol_vm
  all_data['sol_angle']      = sol_angle
  all_data['sol_cvalues']    = sol_cvalues
  all_data['sol_svalues']    = sol_svalues
  all_data['sol_Pfvalues']   = sol_Pfvalues
  all_data['sol_Ptvalues']   = sol_Ptvalues
  all_data['sol_Qfvalues']   = sol_Qfvalues
  all_data['sol_Qtvalues']   = sol_Qtvalues
  all_data['sol_GenPvalues'] = sol_GenPvalues
  all_data['sol_GenQvalues'] = sol_GenQvalues
  
  log.joint(' knitroampl multi-period solution loaded\n')
  

# Writes the current cutting-plane solution to a .txt file
# (variables are sorted by type and index)

def writesol(log,all_data):

  baseMVA      = all_data['baseMVA']
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
  T            = all_data['T']
  alphadic     = all_data['alphadic']
  casetype     = all_data['casetype']
 
  filename = all_data['sols'] + 'CPsol_' + all_data['casename'] + '_' + str(T) + '_' + casetype + '.txt'
  thefile  = open(filename,'w+')
  
  log.joint(' writing solution to ' + filename + '\n')

  # Get machine name and current time
  machinename        = platform.node()
  now                = time.time()

  # Get Gurobi version
  version_info       = gurobipy.gurobi.version()
  version            = '.'.join(map(str, version_info))
  solver_version     = f"Gurobi {version}"

  # System information
  opsystem           = f"{platform.system()} {platform.release()} ({platform.version()})"
  processor          = platform.processor()
  physical_cores     = psutil.cpu_count(logical=False)
  logical_processors = psutil.cpu_count(logical=True)
  cores              = f"{physical_cores} physical cores, {logical_processors} logical processors"
  ram                = f"{round(psutil.virtual_memory().total / (1024 ** 3))} GB RAM"
      
  thefile.write('/CUTPLANEsolution (in p.u.) : ' + all_data['casename'] + '\n')
  thefile.write('/Date : ' + str(time.strftime('%m-%d-%Y %H:%M:%S %Z', time.localtime(now))) + '\n')
  thefile.write('/MachineName : ' + machinename + '\n')
  thefile.write('/Processor : ' + processor + '\n')
  thefile.write('/OS : ' + opsystem + '\n')
  thefile.write('/Cores : ' + cores + '\n')
  thefile.write('/RAM : ' + ram + '\n')
  thefile.write('/Solver : ' + solver_version + '\n')
  thefile.write('/baseMVA : ' + str(baseMVA) + '\n')
  thefile.write('objvalue ' + str(all_data['objval']) + '\n')
  thefile.write('round ' + str(all_data['round']) + '\n')
  thefile.write('time-periods ' + str(all_data['T']) + '\n')
  
  thefile.write('voltages:\n')

  for buscount in buses.keys():
    for k in range(T):
      bus   = buses[buscount]
      v2val = cvalues[k][bus]
      vval  = (v2val)**(0.5)
      line  = 'bus ' + str(buscount) + ' M ' + str(vval) + ' k ' + str(k) + '\n'
    thefile.write(line)

  thefile.write('power flows and cs variables:\n')

  for branchid in branches.keys():
    for k in range(T):
      branch     = branches[branchid]
      f          = branch.f
      t          = branch.t
      count_of_f = IDtoCountmap[f]
      count_of_t = IDtoCountmap[t]
      Pfval      = Pfvalues[k][branch]
      Ptval      = Ptvalues[k][branch]
      Qfval      = Qfvalues[k][branch]
      Qtval      = Qtvalues[k][branch]
      cftval     = cvalues[k][branch]
      sftval     = svalues[k][branch]
      alpha      = alphadic[branch]

      if alpha < all_data['rho_threshold']:
        i2fval     = i2fvalues[k][branch]
        line = 'branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' Pft ' + str(Pfval) + ' Ptf ' + str(Ptval) + ' Qft ' + str(Qfval) + ' Qtf ' + str(Qtval) + ' cft ' + str(cftval) + ' sft ' + str(sftval) + ' i2ft ' + str(i2fval) + ' k ' + str(k) + '\n'
      else:
        line = 'branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' Pft ' + str(Pfval) + ' Ptf ' + str(Ptval) + ' Qft ' + str(Qfval) + ' Qtf ' + str(Qtval) + ' cft ' + str(cftval) + ' sft ' + str(sftval) + ' k ' + str(k) + '\n'
      thefile.write(line)

  thefile.write('generation:\n')

  for genid in gens.keys():
    for k in range(T):
      gen     = gens[genid] 
      nodeID  = gen.nodeID
      GenPval = GenPvalues[k][gen]
      GenQval = GenQvalues[k][gen]
      line_gen = 'genid ' + str(genid) + ' bus ' + str(nodeID) + ' GP ' + str(GenPval) + ' GQ ' + str(GenQval) + ' k ' + str(k) + '\n'
      thefile.write(line_gen)

  thefile.close()

  log.joint(' done writing CUTPLANE solution to .txt file\n\n')

# Writes the current cutting-plane solution to a .sol file
# (variables are not sorted by type)

def writesol_allvars(log,all_data):
  
  baseMVA      = all_data['baseMVA']
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
  T            = all_data['T']
  alphadic     = all_data['alphadic']
  Pd           = all_data['Pd']
  casetype     = all_data['casetype']
    

  filename    = all_data['sols'] + 'CPsol_' + all_data['casename'] + '_' + str(T) + '_' + casetype + '.sol'  
  thefilevars = open(filename,'w+')

  log.joint(' writing solution to ' + filename + '\n')

  # Get machine name and current time
  machinename        = platform.node()
  now                = time.time()

  # Get Gurobi version
  version_info       = gurobipy.gurobi.version()
  version            = '.'.join(map(str, version_info))
  solver_version     = f"Gurobi {version}"

  # System information
  opsystem           = f"{platform.system()} {platform.release()} ({platform.version()})"
  processor          = platform.processor()
  physical_cores     = psutil.cpu_count(logical=False)
  logical_processors = psutil.cpu_count(logical=True)
  cores              = f"{physical_cores} physical cores, {logical_processors} logical processors"
  ram                = f"{round(psutil.virtual_memory().total / (1024 ** 3))} GB RAM"
      
  thefilevars.write('/CUTPLANEsolution (in p.u.) : ' + all_data['casename'] + '\n')
  thefilevars.write('/Date : ' + str(time.strftime('%m-%d-%Y %H:%M:%S %Z', time.localtime(now))) + '\n')
  thefilevars.write('/MachineName : ' + machinename + '\n')
  thefilevars.write('/Processor : ' + processor + '\n')
  thefilevars.write('/OS : ' + opsystem + '\n')
  thefilevars.write('/Cores : ' + cores + '\n')
  thefilevars.write('/RAM : ' + ram + '\n')
  thefilevars.write('/Solver : ' + solver_version + '\n')
  thefilevars.write('/baseMVA : ' + str(baseMVA) + '\n')
  thefilevars.write('/Objvalue ' + str(all_data['objval']) + '\n')
  thefilevars.write('/Time-periods ' + str(all_data['T']) + '\n')

  
  for buscount in buses.keys():
    for k in range(T):
      bus        = buses[buscount]
      f          = bus.nodeID
      v2value    = cvalues[k][bus]        

      v2name     = 'c_' + str(f) + '_' + str(f) + '_' + str(k)
      v2line     = v2name + ' = ' + str(v2value) + '\n'
      thefilevars.write(v2line)

      IPvalue = - Pd[k][bus]
      IQvalue = - bus.Qd
      IPname  = 'IP_' + str(bus.nodeID) + '_' + str(k)
      IQname  = 'IQ_' + str(bus.nodeID) + '_' + str(k)

      for gencounter in bus.genidsbycount:
        if gens[gencounter].status:
          gen      = gens[gencounter]
          IPvalue += GenPvalues[k][gen]
          IQvalue += GenQvalues[k][gen]

      IPline = IPname + ' = ' + str(IPvalue) + '\n'
      thefilevars.write(IPline)

      IQline = IQname + ' = ' + str(IQvalue) + '\n'
      thefilevars.write(IQline)


  for branchid in branches.keys():
    for k in range(T):
      branch     = branches[branchid]
      f          = branch.f
      t          = branch.t
      count_of_f = IDtoCountmap[f]
      count_of_t = IDtoCountmap[t]
      Pfval      = Pfvalues[k][branch]
      Ptval      = Ptvalues[k][branch]
      Qfval      = Qfvalues[k][branch]
      Qtval      = Qtvalues[k][branch]
      cftval     = cvalues[k][branch]
      sftval     = svalues[k][branch]
      alpha      = alphadic[branch]
      
      Pfname  = 'P_' + str(branchid) + '_' + str(f) + '_' + str(t) + '_' + str(k)
      Ptname  = 'P_' + str(branchid) + '_' + str(t) + '_' + str(f) + '_' + str(k)
      Qfname  = 'Q_' + str(branchid) + '_' + str(f) + '_' + str(t) + '_' + str(k)
      Qtname  = 'Q_' + str(branchid) + '_' + str(t) + '_' + str(f) + '_' + str(k)
      cftname = 'c_' + str(branchid) + '_' + str(f) + '_' + str(t) + '_' + str(k)
      sftname = 's_' + str(branchid) + '_' + str(f) + '_' + str(t) + '_' + str(k)

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

      if alpha < all_data['rho_threshold']:
        i2fval     = i2fvalues[k][branch]
        i2fname = 'i2_'+ str(branchid) + '_' + str(f) + '_' + str(t) + '_' + str(k)      
        i2fline  = i2fname + ' = ' + str(i2fval) + '\n'
        thefilevars.write(i2fline)                        


  for genid in gens.keys():
    for k in range(T):
      gen     = gens[genid] 
      nodeID  = gen.nodeID
      GenPval = GenPvalues[k][gen]
      GenQval = GenQvalues[k][gen]
      GPname  = "GP_" + str(genid) + "_" + str(nodeID) + "_" + str(k)
      GQname  = "GQ_" + str(genid) + "_" + str(nodeID) + "_" + str(k)

      GPline  = GPname + ' = ' + str(GenPval) + '\n'
      thefilevars.write(GPline)

      GQline  = GQname + ' = ' + str(GenQval) + '\n'
      thefilevars.write(GQline)


  log.joint(' done writing CUTPLANE allvars solution to .sol file\n\n')

  thefilevars.close()

# Gets the multi-period active power loads from a .txt file

def getloads(log,all_data,loads):

  lines = loads.readlines()
  lenlines = len(lines) - 1 # END                                                                                        
  
  buses        = all_data['buses']
  IDtoCountmap = all_data['IDtoCountmap']
  T            = all_data['T']

  Pd      = {}
  for k in range(T):
    Pd[k] = {}
    
  linenum = 0

  log.joint(' reading file with loads\n')
  while linenum < lenlines:
    thisline = lines[linenum].split()
    buscount         = int(thisline[1])
    bus              = buses[buscount]
    k                = int(thisline[5])
    load             = float(thisline[7])
    Pd[k][bus]       = load
    linenum         += 1

  return  Pd


# Gets the multi-period ramping rates from a .txt file

def getrampr(log,all_data,loads):

  lines = loads.readlines()
  lenlines = len(lines) - 1 # END                                                                                        
  
  gens         = all_data['gens']
  IDtoCountmap = all_data['IDtoCountmap']
  T            = all_data['T']

  rampru      = {}
  ramprd      = {}
  for k in range(T):
    rampru[k] = {}
    ramprd[k] = {}    
    
  linenum = 0

  log.joint(' reading file with ramprates\n')
  while linenum < lenlines:
    thisline = lines[linenum].split()
    gencount         = int(thisline[1])
    gen              = gens[gencount]
    k                = int(thisline[5])
    rpru             = float(thisline[7])
    rprd             = float(thisline[9])
    rampru[k][gen]   = rpru
    ramprd[k][gen]   = rprd
    linenum         += 1

  return  rampru, ramprd

# Retrieves the dual variables corresponding to the 
# active power balance constraints, and compute the difference wrt
# the duals of the previous round

def getduals(log,all_data):

  buses    = all_data['buses']
  themodel = all_data['themodel']
  T        = all_data['T']
  rnd      = all_data['round']
  duals    = {}
  
  for k in range(T):
    for bus in buses.values():
      nodeID     = bus.nodeID
      constrname = "PBaldef"+str(bus.nodeID)+"_"+str(k)
      constr     = themodel.getConstrByName(constrname)
      dual       = constr.Pi
      duals[constrname] = dual 
  
  all_data['duals'][rnd] = duals
  
  if rnd > 1:
    sqdiff = 0
    duals_prev = all_data['duals'][rnd-1]
    for constrname in duals.keys():
      sqdiff += (duals_prev[constrname] - duals[constrname])**2

    dual_diff = round(math.sqrt(sqdiff),4)
    all_data['dual_diff'][rnd] = dual_diff
    log.joint(' dual diff at round  ' + str(rnd) + ' = ' + str(dual_diff) + '\n')

# Prints the difference between consecutive vectors of duals variables
# associated to the active power balance constraints

def print_duals(log,all_data):

  rnd = all_data['round']
  
  log.joint('dual diffs:\n')
  for t in range(2,rnd+1):
    log.joint('(' + str(t) + ',' + str(all_data['dual_diff'][t]) + ')')

  log.joint('\n')

