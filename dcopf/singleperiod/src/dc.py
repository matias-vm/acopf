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

globalalldata = {}


def lpformulator_dc(alldata):
    """Formulate DCOPF model and solve it"""

    log = alldata['log']
    log.joint("\nDC formulation.\n")

    starttime = time.time()

    global globalalldata 
    globalalldata = alldata

    # Create model
    with gp.Env() as env, gp.Model('Single-period DCOPF', env=env) as model:
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

    varcount     = 0
    thetavar     = {}
    Pinjvar      = {}
    Pvar_f       = {}
    GenPvar      = {}


    for j in range(1,numbuses+1):
        bus       = buses[j]
        ubound    = 2*math.pi
        lbound    = -ubound

        thetavar[bus]  = model.addVar(obj = 0.0, lb = lbound, ub = ubound,
                                             name = "theta_"+str(bus.nodeID))
        bus.thetavarind = varcount
        varcount += 1
        
        Plbound = Qlbound = -GRB.INFINITY
        Pubound = Qubound = GRB.INFINITY

        Pubound, Plbound = computebalbounds(log, alldata, bus)

        # Pinjvar is the variable modeling total active power injected by bus j
        # into the branches incident with j
            
        Pinjvar[bus] = model.addVar(obj = 0.0, lb = Plbound, ub = Pubound,
                                           name = "IP_%d"% bus.nodeID )
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
            GenPvar[gen] = model.addVar(obj = 0.0, lb = lower, ub = upper,
                                           name = "GP_%d_%d"%(gen.count,gen.nodeID))
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


        Pvar_f[branch] = model.addVar(obj = 0.0, lb = lbound, ub = ubound,
                                         name = "P_%d_%d_%d"%(j, busf.nodeID, bust.nodeID))
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

        lincostexpr += lincoeff * GenPvar[gen]
        varobjcount +=1

        if gen.costdegree == 2 and gen.costvector[0] != 0:
            qcostexpr   += gen.costvector[0]*GenPvar[gen]*GenPvar[gen]
            varobjcount += 1

    # Constant term
    constexpr   = gp.LinExpr()
    constobjval = 0

    for gen in gens.values():
        if gen.status > 0:
            constobjval += gen.costvector[gen.costdegree]

    constvar   = model.addVar(lb = 1.0, ub = 1.0,name = "constant")
    constexpr += constobjval * constvar

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

        branch.Pfcname = "Pdef_%d_%d_%d"%(j, f, t)

        if branch.status:
            coeff = 1/(branch.x*branch.ratio)
            model.addConstr(Pvar_f[branch] == coeff*thetavar[busf] - coeff*thetavar[bust] - coeff*branch.angle_rad, name = branch.Pfcname)

        else:
            model.addConstr(Pvar_f[branch] == 0, name = branch.Pfcname)

        count += 1
            
    log.joint("    %d active power flow definitions added.\n"%count)
                  

    # Balance constraints
    log.joint("  adding constraints stating bus injection = total outgoing power flow.\n")
    count = 0
    for j in range(1,1+numbuses):
        bus  = buses[j]

        expr = gp.LinExpr()
        for branchid in bus.frombranchids.values():
            expr.add(Pvar_f[branches[branchid]])

        for branchid in bus.tobranchids.values():
            expr.add(-Pvar_f[branches[branchid]])

        model.addConstr(expr == Pinjvar[bus], name = "PBaldef_%d_%d"%(j, bus.nodeID))            
        
        count += 1
    log.joint("    %d constraints added.\n"%count)

    # Injection defs
    log.joint("  adding injection definition constraints.\n")
    count = 0
    for j in range(1,1+numbuses):
        bus  = buses[j]

        expr = gp.LinExpr()
        if len(bus.genidsbycount) > 0:
            for genid in bus.genidsbycount:
                gen = gens[genid]
                expr.add(GenPvar[gen])

        model.addConstr(Pinjvar[bus] == expr - bus.Pd, name = "Bus_PInj_%d_"% j )

        count += 1

    log.joint("    %d injection definition constraints added.n"%count)



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


def computebalbounds(log, alldata, bus):

  # max/min generations                                                                                 
  gens    = alldata['gens']

  Pubound = Plbound = 0
  Qubound = Qlbound = 0

  for gencounter in bus.genidsbycount:
    if gens[gencounter].status:
      Pubound += gens[gencounter].Pmax
      Plbound += gens[gencounter].Pmin

  Pubound -= bus.Pd
  Plbound -= bus.Pd

  if bus.nodetype == 4:
      Pubound = Plbound = 0

  return Pubound, Plbound




