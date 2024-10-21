################################################################################                            
##                                                                            ##                           
## This code was adapted and extended from the Gurobi OptiMod                 ##
## 'Optimal Power Flow', written by Daniel Bienstock. The extension was       ##
## written and it is being maintained by Matias Villagra                      ##
##                                                                            ##                         
## Please report any bugs or issues to                                        ##                           
##                        mjv2153@columbia.edu                                ##                           
##                                                                            ##                         
## Jul 2024                                                                   ##                      
################################################################################


import math
import time
import gurobipy as gp
from gurobipy import GRB
from myutils import break_exit
import numpy as np

globalalldata = {}


def lpformulator_dc_mtp(alldata):
    """Formulate DCOPF model and solve it"""

    log = alldata['log']
    log.joint("\nDC formulation.\n")

    starttime = time.time()

    global globalalldata 
    globalalldata = alldata

    # Create model
    with gp.Env() as env, gp.Model('Multi-period DCOPF', env=env) as model:
        # Add model variables and constraints
        lpformulator_dc_body(alldata, model)

        # Optimize model
        lpformulator_dc_opt(alldata, model)

        endtime = time.time()
        log.joint("Overall time taken (model construction + optimization): %f s.\n"%(endtime-starttime))


def lpformulator_dc_body(alldata, model):
    """Helper function for adding variables and constraints to the model"""

    log           = alldata['log']
    casename      = alldata['casename']
    T             = alldata['T']
    loadsfilename = alldata['loadsfilename']
    rampfilename  = alldata['rampfilename']

    # Getting loads
    loads = open(loadsfilename,"r")
    alldata['Pd'] = getloads(log,alldata,loads)
    log.joint(" Loads obtained\n")

    # Getting ramping rates
    rampr = open(rampfilename,"r")
    rampru, ramprd       = getrampr(log,alldata,rampr)
    log.joint(" Ramp rates obtained\n")
    alldata['rampru'] = rampru
    alldata['ramprd'] = ramprd

    # Create model variables
    lpformulator_dc_create_vars(alldata, model)
    
    # Create model constraints
    lpformulator_dc_create_constraints(alldata, model)

    model.update() 
    
    log.joint("Constructed DCOPF model with %d variables and %d constraints.\n\n"%(model.NumVars, model.NumConstrs))

    if alldata['writeLP']:
        model.write(alldata['lpfilename'])
        log.joint('Wrote LP to ' + alldata['lpfilename'] + '\n')

    alldata['model'] = model

def lpformulator_dc_create_vars(alldata, model):
    """Create model variables for DCOPF"""

    log = alldata['log']
    log.joint('Creating variables.\n')

    numbuses     = alldata['numbuses']
    buses        = alldata['buses']
    IDtoCountmap = alldata['IDtoCountmap']
    gens         = alldata['gens']
    T            = alldata['T']

    varcount     = 0
    thetavar     = {}
    Pinjvar      = {}
    Pvar_f       = {}
    GenPvar      = {}

    for k in range(T):
        thetavar[k] = {}
        Pinjvar[k]  = {}
        Pvar_f[k]   = {}
        GenPvar[k]  = {}


    for j in range(1,numbuses+1):
        bus       = buses[j]
        ubound    = 2*math.pi
        lbound    = -ubound

        
        for k in range(T):
            
            thetavar[k][bus]  = model.addVar(obj = 0.0, lb = lbound, ub = ubound,
                                             name = "theta_"+str(bus.nodeID)+"_"+str(k))
            bus.thetavarind = varcount
            varcount += 1
        
            Plbound = Qlbound = -GRB.INFINITY
            Pubound = Qubound = GRB.INFINITY

            Pubound, Plbound = computebalbounds(log, alldata, bus, k)

            # Pinjvar is the variable modeling total active power injected by bus j
            # into the branches incident with j
            
            Pinjvar[k][bus] = model.addVar(obj = 0.0, lb = Plbound, ub = Pubound,
                                           name = "IP_%d_%d"%(bus.nodeID,k))
            bus.Pinjvarind = varcount
            varcount += 1

            # Generator variables
            for genid in bus.genidsbycount:
                gen   = gens[genid]
                lower = gen.Pmin*gen.status
                upper = gen.Pmax*gen.status
                #if bus.nodetype == 3:
                #  upper = GRB.INFINITY
                #  lower = -GRB.INFINITY  #ignoring slack bus
                GenPvar[k][gen] = model.addVar(obj = 0.0, lb = lower, ub = upper,
                                            name = "GP_%d_%d_%d"%(gen.count,gen.nodeID,k))
                gen.Pvarind = varcount
                varcount += 1
        
    alldata['LP']['thetavar'] = thetavar
    alldata['LP']['Pinjvar']  = Pinjvar
    alldata['LP']['GenPvar']  = GenPvar    

    # Branch related variables
    branches    = alldata['branches']
    numbranches = alldata['numbranches']
    
    for j in range(1,1+numbranches):
        branch     = branches[j]
        f          = branch.f
        t          = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        busf       = buses[count_of_f]
        bust       = buses[count_of_t]
        if branch.constrainedflow:
            ubound     = branch.limit
        else:
            ubound     = alldata['sumPd'] #DC
        lbound     = -ubound

        for k in range(T):
            Pvar_f[k][branch] = model.addVar(obj = 0.0, lb = lbound, ub = ubound,
                                             name = "P_%d_%d_%d_%d"%(j, busf.nodeID, bust.nodeID,k))
            branch.Pftvarind = varcount
            varcount += 1

    alldata['LP']['Pvar_f']     = Pvar_f
    
    # OBJECTIVE

    # Linear and Quadratic terms
    varobjcount = 0

    lincostexpr = gp.LinExpr()
    qcostexpr   = gp.QuadExpr()

    for gen in gens.values():
        lincoeff = gen.costvector[gen.costdegree-1]

        for k in range(T):
            lincostexpr += lincoeff * GenPvar[k][gen]
            varobjcount +=1

            if gen.costdegree == 2 and gen.costvector[0] != 0:
                qcostexpr += gen.costvector[0]*GenPvar[k][gen]*GenPvar[k][gen]
                varobjcount +=1

    # Constant term
    constexpr   = gp.LinExpr()
    constobjval = 0

    for gen in gens.values():
        if gen.status > 0:
            constobjval += gen.costvector[gen.costdegree]

    constvar   = model.addVar(lb = 1.0, ub = 1.0,name = "constant")
    constexpr += T * constobjval * constvar

    varobjcount += 1

    log.joint('  %d terms in the objective\n' %varobjcount)
    
    model.setObjective(constexpr + lincostexpr + qcostexpr)
  
    model.update()

    
def lpformulator_dc_create_constraints(alldata, model):
    """"Create constraint for ACOPF"""

    numbuses     = alldata['numbuses']
    buses        = alldata['buses']
    numbranches  = alldata['numbranches']
    branches     = alldata['branches']
    gens         = alldata['gens']
    IDtoCountmap = alldata['IDtoCountmap']
    log          = alldata['log']
    T            = alldata['T']
    rampru       = alldata['rampru']
    ramprd       = alldata['ramprd']
    Pd           = alldata['Pd']

    thetavar     = alldata['LP']['thetavar']
    Pvar_f       = alldata['LP']['Pvar_f']
    Pinjvar      = alldata['LP']['Pinjvar']
    GenPvar      = alldata['LP']['GenPvar']
    log          = alldata['log']


    log.joint("Creating constraints.\n")

    # Active PF defs
    log.joint("  adding active power flow definitions.\n")
    count = 0
    for j in range(1,1+numbranches):
        branch     = branches[j]
        f          = branch.f
        t          = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        busf       = buses[count_of_f]
        bust       = buses[count_of_t]
        for k in range(T):
            branch.Pfcname = "Pdef_%d_%d_%d_%d"%(j, f, t,k)
            if branch.status:
                coeff = 1/(branch.x*branch.ratio)
                model.addConstr(Pvar_f[k][branch] == coeff*thetavar[k][busf] - coeff*thetavar[k][bust] - coeff*branch.angle_rad, name = branch.Pfcname)

            else:
                model.addConstr(Pvar_f[k][branch] == 0, name = branch.Pfcname)
            count += 1
            
    log.joint("    %d active power flow definitions added.\n"%count)
                  

    # Balance constraints
    log.joint("  adding constraints stating bus injection = total outgoing power flow.\n")
    count = 0
    for j in range(1,1+numbuses):
        bus  = buses[j]
        for k in range(T):
            expr = gp.LinExpr()
            for branchid in bus.frombranchids.values():
                expr.add(Pvar_f[k][branches[branchid]])

            for branchid in bus.tobranchids.values():
                expr.add(-Pvar_f[k][branches[branchid]])

            model.addConstr(expr == Pinjvar[k][bus], name = "PBaldef%d_%d_%d"%(j, bus.nodeID,k))            
        
            count += 1
    log.joint("    %d constraints added.\n"%count)

    # Injection defs
    log.joint("  adding injection definition constraints.\n")
    count = 0
    for j in range(1,1+numbuses):
        bus  = buses[j]
        for k in range(T):
            expr = gp.LinExpr()
            if len(bus.genidsbycount) > 0:
                for genid in bus.genidsbycount:
                    gen = gens[genid]
                    expr.add(GenPvar[k][gen])

            model.addConstr(Pinjvar[k][bus] == expr - Pd[k][bus], name = "Bus_PInj_%d_%d"%(j,k))
            count += 1

    log.joint("    %d injection definition constraints added.n"%count)

    # Ramping contraints                                                                                                  
    log.joint('  adding ramping constraints...\n')

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

            abs_gen[k][gen] = model.addVar(lb = 0, name = 'abs_gen_'+str(genid)+'_'+str(k))
            model.addConstr(GenPvar[k+1][gen] - GenPvar[k][gen] - rpu * abs_gen[k][gen] <= 0, name = constrname_rup)
            model.addConstr(GenPvar[k][gen] - rpd * abs_gen[k][gen] - GenPvar[k+1][gen] <= 0, name = constrname_rdown)

            model.addConstr(GenPvar[k][gen] - abs_gen[k][gen] <= 0, name = constrname_rup + '_1')
            model.addConstr(- GenPvar[k][gen] - abs_gen[k][gen] <= 0, name = constrname_rup + '_2')

            count += 4


def lpformulator_dc_opt(alldata, model):

    log          = alldata['log']
    numbranches  = alldata['numbranches']
    branches     = alldata['branches']

    model.Params.OptimalityTol  = 1.0e-8
    model.Params.FeasibilityTol = 1.0e-8 
    model.Params.NumericFocus   = 1
    model.Params.BarHomogeneous = 1
    model.Params.method    = 2
    model.Params.Crossover = 0
    feastol = model.Params.FeasibilityTol
    opttol  = model.Params.OptimalityTol

    model.optimize()
        
    # Check model status and re-optimize or try computing an IIS if necessary
    if model.status == GRB.INF_OR_UNBD:
        log.joint("\nModel Status: infeasible or unbounded.\n")
        log.joint("Re-optimizing with DualReductions turned off.\n\n")
        model.Params.DualReductions = 0
        model.optimize()
        
    if model.status == GRB.INFEASIBLE:
        log.joint("\nModel Status: infeasible.")
        log.joint("Computing IIS...\n\n")
        model.computeIIS()
        log.joint("\nIIS computed, writing IIS to file dcopf_mtp.ilp\n\n")
        model.write("acopfmodel.ilp")
        
    elif model.status == GRB.UNBOUNDED:
        log.joint("\nModel Status: unbounded.\n\n")

    elif model.status == GRB.INTERRUPTED:
        log.joint("\nModel Status: interrupted.\n\n")

    elif model.status == GRB.OPTIMAL:
        log.joint("\nModel Status: optimal.\n\n")

    alldata['optstatus']  = model.status
    alldata['numvars']    = model.numvars
    alldata['numconstrs'] = model.numconstrs
    
    if model.SolCount > 0:
            log.joint("Objective value = %g"%model.objVal)
            model.printQuality()
            alldata['objval']     = model.objVal
            if model.status == 2 or model.status == 13:
                alldata['dinfs_scaled'] = model.DualResidual
                alldata['dinfs'] = model.DualSResidual
            else:
                alldata['dinfs_scaled'] = -1
                alldata['dinfs']        = -1
    else:
        alldata['objval']     = -1
        alldata['dinfs_scaled'] = -1
        alldata['dinfs']        = -1

    try:
        alldata['Kappa']      = model.Kappa
        alldata['KappaExact'] = model.KappaExact
    except:
        alldata['Kappa']      = -1
        alldata['KappaExact'] = -1


def computebalbounds(log, alldata, bus, k):

  # max/min generations                                                                                 
  gens    = alldata['gens']
  Pd      = alldata['Pd']

  Pubound = Plbound = 0
  Qubound = Qlbound = 0

  for gencounter in bus.genidsbycount:
    if gens[gencounter].status:
      Pubound += gens[gencounter].Pmax
      Plbound += gens[gencounter].Pmin


  Pubound -= Pd[k][bus]
  Plbound -= Pd[k][bus]

  if bus.nodetype == 4:
      Pubound = Plbound = 0
      Pd[k][bus] = 0

  return Pubound, Plbound

def getloads(log,alldata,loads):

  lines = loads.readlines()
  lenlines = len(lines) - 1 # END                                                                        
  buses        = alldata['buses']
  IDtoCountmap = alldata['IDtoCountmap']
  T            = alldata['T']

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


def getrampr(log,alldata,loads):

  lines = loads.readlines()
  lenlines = len(lines) - 1 # END                                                                         
  gens         = alldata['gens']
  IDtoCountmap = alldata['IDtoCountmap']
  T            = alldata['T']

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
