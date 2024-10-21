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

from amplpy import AMPL
from myutils import breakexit
from log import danoLogger
import time
import math
import os
import platform
import psutil

def goac_mtp2(log,all_data):
    
    log.joint(' creating ampl object ...\n')

    ampl         = AMPL()
    casename     = all_data['casename']
    casetype     = all_data['casetype']
    modfile      = all_data['modfile']
    solver       = all_data['solver']
    IDtoCountmap = all_data['IDtoCountmap']
    T            = all_data['T']
    
    log.joint(" reading modfile ...\n")
    
    t0 = time.time()
    ampl.read('../modfiles/'+modfile)
    t1 = time.time()
    log.joint(" modfile read in time " + str(t1-t0))

    # Setting up AMPL options
    ampl.setOption('display_precision',0)
    ampl.setOption('expand_precision',0)
    ampl.setOption('show_stats',2)
    ampl.setOption('solver',solver)    

    log.joint("\n solver set to " + solver + "\n")
    
    ampl.setOption('presolve',0)        
    log.joint(' AMPL presolve off\n')
        
    # Solver options
    feastol_abs          = all_data['feastol_abs']
    opttol_abs           = all_data['opttol_abs']
    feastol_rel          = all_data['feastol_rel']
    opttol_rel           = all_data['opttol_rel']
    honorbnds            = all_data['honorbnds']
    ftol                 = all_data['ftol']
    scale                = all_data['scale']
    ftol_iters           = all_data['ftol_iters']
    linsolver            = all_data['linsolver']
    max_time             = all_data['max_time']
    wstart               = all_data['wstart']
    bar_initmu           = all_data['bar_initmu']
    blasoption           = all_data['blasoption']
    linsolver_numthreads = all_data['linsolver_numthreads']
    blas_numthreads      = all_data['blas_numthreads']
    bar_murule           = all_data['bar_murule']
    mylogfile            = all_data['mylogfile']
    
    solver_options = "option knitro_options 'feastol=" + feastol_rel + " feastol_abs=" + feastol_abs + " opttol=" + opttol_rel + " opttol_abs=" + opttol_abs + " ftol=" + ftol + " ftol_iters=" + ftol_iters + " honorbnds=" + honorbnds + " blasoption=" + blasoption + " blas_numthreads=" + blas_numthreads + " linsolver_numthreads=" + linsolver_numthreads + " numthreads=20 linsolver=" + linsolver + " maxtime_real=" + max_time + " strat_warm_start=" + wstart + " bar_murule=" + bar_murule + " bar_initmu=" + bar_initmu + " outname=" + mylogfile + " outmode=2';"
    
    ampl.eval(solver_options)

    if all_data['multistart']:
        solver_options = "option knitro_options 'feastol_abs=" + feastol + " opttol_abs=" + opttol + " blasoptionlib=1 numthreads=40 linsolver=" + linsolver + " maxtime_real=" + max_time + " strat_warm_start=" + wstart + " bar_initmu=" + bar_initmu + " ms_enable=1 ms_numthreads=20 ms_maxsolves=4 ms_terminate =1';"
        ampl.eval(solver_options)
        
    # Initial solution (flat-start or previously computed AC solutions)
    Vinit_mtp, thetadiffinit_mtp = initial_solution_mtp2(log,all_data)

    # Getting multi-time period loads and ramp rates from file
    Pd_mtp                 = getloads(log,all_data)
    rampru_mtp, ramprd_mtp = getrampr(log,all_data)
    all_data['Pd_mtp']     = Pd_mtp

    # Initializing Pmax_mtp and Pmin_mtp
    Pmax_mtp, Pmin_mtp = datastruct_Pmax_min(log,all_data)
    
    # Setting up buses
    buses        = {}
    Pd           = {}
    Qd           = {}
    Vmax         = {}
    Vmin         = {}    
    branches_f   = {}
    branches_t   = {}
    bus_gens     = {}
    bus_Gs       = {}
    bus_Bs       = {}
    theta_min    = {}
    theta_max    = {}
    Vinit        = {}

    
    for bus in all_data['buses'].values():
        buscount = bus.count

        if ( len(bus.frombranchids) == 0 ) and ( len(bus.tobranchids) == 0 ): 
            log.joint(' isolated bus ' + str(bus.nodeID) + ' busid ' + str(buscount) + '\n') 
        
        buses[buscount]      = bus.nodeID 
        Qd[buscount]         = bus.Qd
        Vmin[buscount]       = bus.Vmin
        Vmax[buscount]       = bus.Vmax
        branches_f[buscount] = []
        branches_t[buscount] = []
        bus_gens[buscount]   = bus.genidsbycount 

        if buscount == all_data['refbus']:
            theta_min[buscount] = 0
            theta_max[buscount] = 0
        else:
            theta_min[buscount] = -2 * math.pi
            theta_max[buscount] = 2 * math.pi

        if ( bus.Gs != 0 ) and ( ( len(bus.frombranchids) != 0 ) or ( len(bus.tobranchids) != 0 ) ):
            bus_Gs[buscount] = bus.Gs
        if ( bus.Bs != 0 ) and ( ( len(bus.frombranchids) != 0 ) or ( len(bus.tobranchids) != 0 ) ):
            bus_Bs[buscount] = bus.Bs                

        for branchid in bus.frombranchids.values():
            branches_f[buscount].append(branchid) 
            
        for branchid in bus.tobranchids.values():
            branches_t[buscount].append(branchid)      

            
    ampl.getSet('buses').setValues(list(buses))
    ampl.getSet('bus_Gs').setValues(list(bus_Gs))
    ampl.getSet('bus_Bs').setValues(list(bus_Bs))
    ampl.get_parameter('Gs').setValues(list(bus_Gs.values()))
    ampl.get_parameter('Bs').setValues(list(bus_Bs.values()))
    ampl.get_parameter('Qd').set_values(Qd)
    ampl.get_parameter('Vmax').set_values(Vmax)
    ampl.get_parameter('Vmin').set_values(Vmin)
    ampl.get_parameter('theta_min').setValues(theta_min)
    ampl.get_parameter('theta_max').setValues(theta_max)

    branches_f_set = ampl.getSet('branches_f')
    branches_t_set = ampl.getSet('branches_t')
    bus_gens_set   = ampl.getSet('bus_gens')
    
    for bus in all_data['buses'].values():
        branches_f_set[bus.count].setValues(branches_f[bus.count])
        branches_t_set[bus.count].setValues(branches_t[bus.count])
        bus_gens_set[bus.count].setValues(bus_gens[bus.count])
    
    # Setting up branches
    branches      = {}
    U             = {}
    Gtt           = {}
    Btt           = {}
    Gff           = {}
    Bff           = {}
    Gtf           = {}
    Btf           = {}
    Gft           = {}
    Bft           = {}
    bus_f         = {}
    bus_t         = {}
    maxangle      = {}
    minangle      = {}
    thetadiffinit = {}
    
    for branch in all_data['branches'].values():
        branchcount            = branch.count

        if branch.status != 1: 
            log.joint(' branch ' + str(branchcount) + ' OFF, we skip it\n')
            breakexit('check')
            continue

        branches[branchcount]  = (branch.id_f,branch.id_t)
        U[branchcount]         = branch.limit
        Gtt[branchcount]       = branch.Gtt
        Btt[branchcount]       = branch.Btt
        Gff[branchcount]       = branch.Gff
        Bff[branchcount]       = branch.Bff
        Gtf[branchcount]       = branch.Gtf
        Btf[branchcount]       = branch.Btf
        Gft[branchcount]       = branch.Gft
        Bft[branchcount]       = branch.Bft
        bus_f[branchcount]     = branch.id_f
        bus_t[branchcount]     = branch.id_t
        maxangle[branchcount]  = branch.maxangle_rad
        minangle[branchcount]  = branch.minangle_rad
        
            
    ampl.getSet('branches').setValues(list(branches))
    ampl.get_parameter('Gtt').setValues(Gtt)
    ampl.get_parameter('Btt').setValues(Btt)
    ampl.get_parameter('Gff').setValues(Gff)
    ampl.get_parameter('Bff').setValues(Bff)
    ampl.get_parameter('Gtf').setValues(Gtf)
    ampl.get_parameter('Btf').setValues(Btf)
    ampl.get_parameter('Gft').setValues(Gft)
    ampl.get_parameter('Bft').setValues(Bft) 
    ampl.get_parameter('U').setValues(U)
    ampl.get_parameter('bus_f').setValues(bus_f)
    ampl.get_parameter('bus_t').setValues(bus_t)
    ampl.get_parameter('maxangle').setValues(maxangle)
    ampl.get_parameter('minangle').setValues(minangle)


    # Setting up generators
    gens      = {}
    Pmax      = {}
    Pmin      = {}
    Qmax      = {}
    Qmin      = {}
    fixedcost = {}
    lincost   = {}
    quadcost  = {} 

    for gen in all_data['gens'].values():
        gencount            = gen.count 
        gens[gencount]      = gen.nodeID
        Qmax[gencount]      = gen.Qmax * gen.status
        Qmin[gencount]      = gen.Qmin * gen.status

        buscount = IDtoCountmap[gen.nodeID]
        bus      = all_data['buses'][buscount]
        if bus.nodetype == 3:
            Qmax[gencount] = 2 * all_data['summaxgenQ']
            Qmin[gencount] = - 2 * all_data['summaxgenQ']
    
        
        fixedcost[gencount] = gen.costvector[2] * gen.status
        lincost[gencount]   = gen.costvector[1]
        quadcost[gencount]  = gen.costvector[0]
        
        
    ampl.getSet('gens').setValues(list(gens))
    ampl.get_parameter('Qmax').setValues(Qmax)
    ampl.get_parameter('Qmin').setValues(Qmin)
    ampl.get_parameter('fixedcost').setValues(fixedcost)
    ampl.get_parameter('lincost').setValues(lincost)
    ampl.get_parameter('quadcost').setValues(quadcost)
    
    log.joint(" sets and time-invariant parameters loaded\n")

    # Initializing data structure for AC solutions
    datastruct_acsolutions(log,all_data)

    # Main loop
    objvalue = 0
    
    for k in range(T):

        # Setting up active-power demand and initial voltages
        for buscount in all_data['buses'].keys():
            Pd[buscount]    = Pd_mtp[k][buscount]
            Vinit[buscount] = Vinit_mtp[k][buscount]

        ampl.get_parameter('Pd').set_values(Pd)
        ampl.get_parameter('Vinit').set_values(Vinit)

        
        # Setting up initial angel differences
        for branchcount in all_data['branches'].keys():
            thetadiffinit[branchcount] = thetadiffinit_mtp[k][branchcount]
            
        ampl.get_parameter('thetadiffinit').set_values(thetadiffinit)

        
        # Setting up max/min active power generation
        for gen in all_data['gens'].values():
            gencount = gen.count
            nodeID  = gen.nodeID            
            if k == 0:
                Pmax[gencount] = gen.Pmax * gen.status
                Pmin[gencount] = gen.Pmin * gen.status
            else:
                Pmax[gencount] = Pmax_mtp[k-1][gencount] * gen.status
                Pmin[gencount] = Pmin_mtp[k-1][gencount] * gen.status

        ampl.get_parameter('Pmax').setValues(Pmax)
        ampl.get_parameter('Pmin').setValues(Pmin)

        
        # Expand
        if all_data['expand']:
            filename = casename + "_" + casetype + "_" + str(k) + "_acNLP.out"
            log.joint('Now expanding to %s.\n'%(filename))
            amplstate = 'expand; display {j in 1.._nvars} (_varname[j],_var[j].lb,_var[j].ub);' 
            modelout = ampl.getOutput(amplstate)
            outfile = open(filename,"w")
            outfile.write("model = " + str(modelout) + "\n")
            outfile.close()

        # Reset variables for flat start
        amplstate='reset data v; reset data thetadiff;'
        ampl.eval(amplstate)        
        
        amplstate='for {i in buses}{ let v[i] := 1; };'
        ampl.eval(amplstate)
        amplstate='for {i in branches}{ let thetadiff[i] := 0; };'        
        ampl.eval(amplstate)            
            
        # Solve
        log.joint(" solving model ...\n")
        t0 = time.time()
        ampl.solve()
        t1 = time.time()
        
        # Updating Knitro options
        bar_initmu = "1e-2"
        solver_options = "option knitro_options 'feastol=" + feastol_rel + " feastol_abs=" + feastol_abs + " opttol=" + opttol_rel + " opttol_abs=" + opttol_abs + " ftol=" + ftol + " ftol_iters=" + ftol_iters + " honorbnds=" + honorbnds + " blasoption=" + blasoption + " blas_numthreads=" + blas_numthreads + " linsolver_numthreads=" + linsolver_numthreads + " numthreads=20 linsolver=" + linsolver + " maxtime_real=" + max_time + " strat_warm_start=" + wstart + " bar_murule=" + bar_murule + " bar_initmu=" + bar_initmu + "';"
        #options   = "option knitro_options 'feastol_abs=" + feastol_abs + " opttol_abs=" + opttol_abs + " blasoptionlib=1 numthreads=20 linsolver=" + linsolver + " maxtime_real=" + max_time + " strat_warm_start=" + wstart + " bar_initmu=" + bar_initmu + "';"
        ampl.eval(solver_options)


        log.joint("\n ===============================================================\n")
        log.joint(" ===============================================================\n")
    
        # Get solution
        solver_status = ampl.getOutput("display solve_result_num;").split(' ')[2].strip()
        
        total_costvar              = ampl.get_objective("total_cost")
        vvar                       = ampl.get_variable("v")
        thetavar                   = ampl.get_variable("theta")
        Pfvar                      = ampl.get_variable("Pf")
        Ptvar                      = ampl.get_variable("Pt")
        Qfvar                      = ampl.get_variable("Qf")
        Qtvar                      = ampl.get_variable("Qt")
        GenPvar                    = ampl.get_variable("Pg")
        GenQvar                    = ampl.get_variable("Qg")
        all_data['objvalue'][k]    = total_costvar.get().value()
        all_data['vvalues'][k]     = vvar.get_values().to_dict()
        all_data['thetavalues'][k] = thetavar.get_values().to_dict()
        all_data['GenPvalues'][k]  = GenPvar.get_values().to_dict()
        all_data['GenQvalues'][k]  = GenQvar.get_values().to_dict()
        all_data['Pfvalues'][k]    = Pfvar.get_values().to_dict()
        all_data['Ptvalues'][k]    = Ptvar.get_values().to_dict()
        all_data['Qfvalues'][k]    = Qfvar.get_values().to_dict()
        all_data['Qtvalues'][k]    = Qtvar.get_values().to_dict()
           
        sumPd = 0
        sumQd = 0
        for buscount in Pd.keys():
            sumPd += Pd[buscount]
            sumQd += Qd[buscount]

        PLoss    = {}
        QLoss    = {}
        sumPLoss = 0
        sumQLoss = 0

        for branchcount in all_data['Pfvalues'][k].keys():
            PLoss[branchcount] = all_data['Pfvalues'][k][branchcount] + all_data['Ptvalues'][k][branchcount]
            QLoss[branchcount] = all_data['Qfvalues'][k][branchcount] + all_data['Qtvalues'][k][branchcount]
            sumPLoss += PLoss[branchcount]
            sumQLoss += QLoss[branchcount]

        sumGenP = 0
        sumGenQ = 0
        for gencount in all_data['GenPvalues'][k].keys():
            sumGenP += all_data['GenPvalues'][k][gencount]
            sumGenP += all_data['GenQvalues'][k][gencount]

        
        all_data['sumPd'][k]         = sumPd
        all_data['sumQd'][k]         = sumQd
        all_data['sumPLoss'][k]      = sumPLoss
        all_data['sumQLoss'][k]      = sumQLoss
        all_data['sumGenP'][k]       = sumGenP
        all_data['sumGenQ'][k]       = sumGenQ      
        all_data['solver_status'][k] = int(solver_status)
        
        timesofar = time.time() - all_data['T0']        

        # Logging solution
        log.joint(' case ' + casename + ' period ' + str(k) + '/'
                  + str(T) + ' ' + casetype + '\n')
        log.joint(' solver ' + all_data['solver'] + ' status '
                  + str(all_data['solver_status'][k]) + '\n')
        log.joint(' objective ' + str(all_data['objvalue'][k]) + '\n')        
        log.joint(" modfile " + all_data['modfile'] + "\n")
        log.joint(" active power generation " + str(sumGenP) + '\n')
        log.joint(" active power demand " + str(sumPd) + '\n')
        log.joint(" active power loss " + str(sumPLoss) + '\n')
        log.joint(" reactive power generation " + str(sumGenQ) + '\n')
        log.joint(" reactive power demand " + str(sumQd) + '\n')
        log.joint(" reactive power loss " + str(sumQLoss) + '\n')
        log.joint(" solver runtime " + str(t1-t0) + '\n')
        log.joint(" time so far " + str(timesofar) + '\n')


        summary_ac = open("summary_ac.log","a+")
        
        #if all_data['solver_status'][k] != 0:
        #    log.joint(' No feasible solution, H2 unsuccesful\n')
        #    break
        #else:
            
        objvalue += all_data['objvalue'][k]
        log.joint(" writing casename, modfile, obj and runtime to summary_ac.log\n")
        summary_ac.write(' case ' + all_data['casename'] + ' casetype ' + all_data['casetype']
                         + ' modfile ' + all_data['modfile'] + ' solver ' + all_data['solver']
                         + ' solver_status ' + str(all_data['solver_status'].values())
                         + ' obj ' + str(all_data['objvalue'][k]) + ' runtime ' + str(timesofar)
                         + '\n')
        summary_ac.close()


        if (all_data['solver_status'][k] != 0 and all_data['solver_status'][k] != 103
            and all_data['solver_status'][k] != 102 and all_data['solver_status'][k] != 101):
            log.joint(' Solution is not feasible, H2 failed\n')
            all_data['h2']['objvalue'] = round(objvalue,2)
            all_data['h2']['time']     = round(timesofar,2)
            all_data['h2']['status']   = all_data['solver_status'][k]
            return 1
        
        log.joint(" ===============================================================\n")
        log.joint(" ===============================================================\n")

        
        # Updating Pmax and Pmin for next time-period
        log.joint(' Updating Pmax and Pmin for next time-period\n')

        if k < T-1:
            for genid in all_data['gens'].keys():
                gen       = all_data['gens'][genid]
                nodeID    = gen.nodeID
                genval    = all_data['GenPvalues'][k][genid]
                absgenval = math.fabs(genval)
                pmax      = gen.status * gen.Pmax
                pmin      = gen.status * gen.Pmin
                Pmax_mtp[k][genid] = min(genval + rampru_mtp[k][genid] * absgenval,pmax)
                Pmin_mtp[k][genid] = max(genval - ramprd_mtp[k][genid] * absgenval,pmin)
                #if genval < 0:
                #    log.joint(' gen ' + str(genid) + ' at bus ' + str(nodeID)
                #              + ' GP ' + str(genval) + '\n')
        
        if all_data['writesol']:
            writesol_k(log,all_data,k)
            writesol_qcqp_allvars_k(log,all_data,k)

    log.joint("\n\n ===============================================================\n")
    log.joint(" ===============================================================\n")
        
    sumPd    = sum(all_data['sumPd'].values())
    sumQd    = sum(all_data['sumQd'].values())    
    sumPLoss = sum(all_data['sumPLoss'].values())
    sumQLoss = sum(all_data['sumQLoss'].values())
    sumGenP  = sum(all_data['sumGenP'].values())
    sumGenQ  = sum(all_data['sumGenQ'].values())    
    status   = sum(all_data['solver_status'].values())
    
    timesofar = time.time() - all_data['T0']

    if status != 0:
        status   = ""
        for k in range(len(all_data['solver_status'])):
            if k < len(all_data['solver_status']) - 1:
                status += str(all_data['solver_status'][k]) + "-"
            else:
                status += str(all_data['solver_status'][k])
    else:
        status = "0"
        
    log.joint(' case ' + casename + ' sum of ' + str(T) + ' periods '
              + casetype + '\n')
    if status == 0:
        log.joint(" solver " + all_data['solver']
                  + " primal-bound found"
                  + " warmstart "
                  + str(all_data['wstart']) + "\n")
    else:
        log.joint(" solver " + all_data['solver']
                  + " warmstart "
                  + str(all_data['wstart'])
                  + " list of error codes: " + status + "\n")
    log.joint(" objective " + str(objvalue) + '\n')        
    log.joint(" modfile " + all_data['modfile'] + "\n")
    log.joint(" active power generation " + str(sumGenP) + '\n')
    log.joint(" active power demand " + str(sumPd) + '\n')
    log.joint(" active power loss " + str(sumPLoss) + '\n')
    log.joint(" reactive power generation " + str(sumGenQ) + '\n')
    log.joint(" reactive power demand " + str(sumQd) + '\n')
    log.joint(" reactive power loss " + str(sumQLoss) + '\n')
    log.joint(" time so far " + str(timesofar) + '\n')


    all_data['h2']['objvalue'] = round(objvalue,2)
    all_data['h2']['time']     = round(timesofar,2)
    all_data['h2']['status']   = status
    
    log.joint(" ===============================================================\n")
    log.joint(" ===============================================================\n")

    
    log.joint(" writing casename, modfile, obj and runtime to summary_ac.log\n")

    summary_ac = open("summary_ac.log","a+")
    
    summary_ac.write(' case ' + all_data['casename'] + ' T ' + str(all_data['T'])
                     + ' casetype ' + all_data['casetype']
                     + ' modfile ' + all_data['modfile'] + ' solver ' + all_data['solver']
                     + ' solver_status ' + status + ' obj ' + str(objvalue) + ' runtime '
                     + str(timesofar) + '\n')
    
    summary_ac.close()

    if all_data['writesol']:
        writesol(log,all_data)
        writesol_qcqp_allvars(log,all_data)
 
    if status == "0 ":
        return 0
    else:
        return 1
        
def writesol(log,all_data):

    casename     = all_data['casename']
    casetype     = all_data['casetype']
    branches     = all_data['branches']
    buses        = all_data['buses']
    gens         = all_data['gens']
    IDtoCountmap = all_data['IDtoCountmap']
    objvalue     = all_data['objvalue']
    vvalues      = all_data['vvalues']
    thetavalues  = all_data['thetavalues']
    GenPvalues   = all_data['GenPvalues']
    GenQvalues   = all_data['GenQvalues']
    Pfvalues     = all_data['Pfvalues']
    Ptvalues     = all_data['Ptvalues']
    Qfvalues     = all_data['Qfvalues']
    Qtvalues     = all_data['Qtvalues']
    T            = all_data['T']
    casetype     = all_data['casetype']


    filename      = all_data['sols'] + 'ACsol_' + casename + '_' + str(T) + '_' + casetype + '.txt'
    thefile       = open(filename,'w+')

    log.joint('\n writing solution to ' + filename + '\n')

    machinename = platform.node()
    now = time.time()
    AMPL_version = 'Version 20190223'
    solver_version = 'Artelys Knitro 14.1.0'
    opsystem = "{} {} ({})".format(platform.system(), platform.release(), platform.version())
    processor = platform.processor() 
    physical_cores = psutil.cpu_count(logical=False)
    logical_processors = psutil.cpu_count(logical=True)
    cores = f"{physical_cores} physical cores, {logical_processors} logical processors"
    ram = f"{round(psutil.virtual_memory().total / (1024 ** 3))} GB RAM"

    thefile.write('/ACsolution : ' + casename + " T" + str(T) + " " + casetype + '\n')
    thefile.write('/Date : ' + str(time.strftime('%m-%d-%Y %H:%M:%S %Z', time.localtime(now))) + '\n')
    thefile.write('/MachineName : ' + machinename + '\n')
    thefile.write('/Processor : ' + processor + '\n')
    thefile.write('/OS : ' + opsystem + '\n')
    thefile.write('/Cores : ' + cores + '\n')
    thefile.write('/RAM : ' + ram + '\n')
    thefile.write('/AMPL : ' + AMPL_version + '\n')
    thefile.write('/Solver : ' + solver_version + '\n')
    thefile.write('objvalue ' + str(sum(all_data['objvalue'].values())) + '\n')
    thefile.write('time-periods ' + str(T) + '\n')
        
    thefile.write('voltages and angles:\n')

    for buscount in buses.keys():
        for k in range(T):
            vval       = vvalues[k][buscount]
            thetaval   = thetavalues[k][buscount]
            line = 'bus ' + str(buscount) + ' M ' + str(vval) + ' A ' + str(thetaval) + ' k ' + str(k) + '\n'
            thefile.write(line)
            
    thefile.write('power flows and cs variables:\n')
    
    for branchid in branches.keys():
        for k in range(T):
            branch     = branches[branchid]
            f          = branch.f
            t          = branch.t
            count_of_f = IDtoCountmap[f]
            count_of_t = IDtoCountmap[t]
            Pfval      = Pfvalues[k][branchid]
            Ptval      = Ptvalues[k][branchid]
            Qfval      = Qfvalues[k][branchid]
            Qtval      = Qtvalues[k][branchid]
            thetaval   = thetavalues[k][count_of_f] - thetavalues[k][count_of_t]
            cftval     = vvalues[k][count_of_f] * vvalues[k][count_of_t] * math.cos(thetaval)
            sftval     = vvalues[k][count_of_f] * vvalues[k][count_of_t] * math.sin(thetaval)        
        
            line = 'branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' Pft ' + str(Pfval) + ' Ptf ' + str(Ptval) + ' Qft ' + str(Qfval) + ' Qtf ' + str(Qtval) + ' cft ' + str(cftval) + ' sft ' + str(sftval) + ' k ' + str(k) + '\n'
            thefile.write(line)

    thefile.write('generation:\n')
        
    for genid in gens.keys():
        for k in range(T):
            gen     = gens[genid] 
            nodeID  = gen.nodeID
            line_gen = 'genid ' + str(genid) + ' bus ' + str(nodeID) + ' GP ' + str(GenPvalues[k][genid]) + ' GQ ' + str(GenQvalues[k][genid]) + ' k ' + str(k) + '\n'
            thefile.write(line_gen)
        
        
    thefile.close()

    log.joint(' done writing knitro solution to .txt file\n')


def writesol_qcqp_allvars(log,all_data):

    casename      = all_data['casename']
    casetype      = all_data['casetype']    
    branches      = all_data['branches']
    buses         = all_data['buses']
    gens          = all_data['gens']
    IDtoCountmap  = all_data['IDtoCountmap']
    objvalue      = all_data['objvalue']
    vvalues       = all_data['vvalues']
    thetavalues   = all_data['thetavalues']
    GenPvalues    = all_data['GenPvalues']
    GenQvalues    = all_data['GenQvalues']
    Pfvalues      = all_data['Pfvalues']
    Ptvalues      = all_data['Ptvalues']
    Qfvalues      = all_data['Qfvalues']
    Qtvalues      = all_data['Qtvalues']
    Pd            = all_data['Pd_mtp']
    T             = all_data['T']
    
    filenamevars  = all_data['sols'] + 'ACsol_' + all_data['casename'] + '_' + str(T) + '_' + casetype + '.sol'
    thefilevars   = open(filenamevars,'w+')
    
    log.joint('\n writing solution to ' + filenamevars + '\n')

    machinename = platform.node()
    now = time.time()
    AMPL_version = 'Version 20190223'
    solver_version = 'Artelys Knitro 14.1.0'
    opsystem = "{} {} ({})".format(platform.system(), platform.release(), platform.version())
    processor = platform.processor()
    physical_cores = psutil.cpu_count(logical=False)
    logical_processors = psutil.cpu_count(logical=True)
    cores = f"{physical_cores} physical cores, {logical_processors} logical processors"
    ram = f"{round(psutil.virtual_memory().total / (1024 ** 3))} GB RAM"

    thefilevars.write('/ACsolution : ' + casename + " T" + str(T) + " " + casetype + '\n')
    thefilevars.write('/Date : ' + str(time.strftime('%m-%d-%Y %H:%M:%S %Z', time.localtime(now))) + '\n')
    thefilevars.write('/MachineName : ' + machinename + '\n')
    thefilevars.write('/Processor : ' + processor + '\n')
    thefilevars.write('/OS : ' + opsystem + '\n')
    thefilevars.write('/Cores : ' + cores + '\n')
    thefilevars.write('/RAM : ' + ram + '\n')
    thefilevars.write('/AMPL : ' + AMPL_version + '\n')
    thefilevars.write('/Solver : ' + solver_version + '\n')
    thefilevars.write('/Objvalue : ' + str(sum(all_data['objvalue'].values())) + '\n')
    thefilevars.write('/Time-periods ' + str(T) + '\n')
    
    for buscount in buses.keys():
        for k in range(T):
            bus        = buses[buscount]
            vval       = vvalues[k][buscount]
            f          = bus.nodeID
            thetaval   = thetavalues[k][buscount]
            v2value    = vval**2
            evalue     = vval * math.cos(thetaval)
            fvalue     = vval * math.sin(thetaval)

            v2name     = 'c_' + str(f) + '_' + str(f) + '_' + str(k)
            ename      = 'e_' + str(f) + '_' + str(k)
            fname      = 'f_' + str(f) + '_' + str(k)

            v2line     = v2name + ' = ' + str(v2value) + '\n'
            thefilevars.write(v2line)

            eline     = ename + ' = ' + str(evalue) + '\n'
            thefilevars.write(eline)

            fline     = fname + ' = ' + str(fvalue) + '\n'
            thefilevars.write(fline)        

            IPvalue = - Pd[k][buscount]
            IQvalue = - bus.Qd
            IPname  = 'IP_' + str(bus.nodeID) + '_' + str(k)
            IQname  = 'IQ_' + str(bus.nodeID) + '_' + str(k)

            for gencounter in bus.genidsbycount:
                if gens[gencounter].status:
                    IPvalue += GenPvalues[k][gencounter]
                    IQvalue += GenQvalues[k][gencounter]

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
            thetaftval = thetavalues[k][count_of_f] - thetavalues[k][count_of_t]
            Pfval      = Pfvalues[k][branchid]
            Ptval      = Ptvalues[k][branchid]
            Qfval      = Qfvalues[k][branchid]
            Qtval      = Qtvalues[k][branchid]
            cftval     = vvalues[k][count_of_f] * vvalues[k][count_of_t] * math.cos(thetaftval)
            sftval     = vvalues[k][count_of_f] * vvalues[k][count_of_t] * math.sin(thetaftval)
            
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

        
    for genid in gens.keys():
        for k in range(T):
            gen     = gens[genid] 
            nodeID  = gen.nodeID
            GenPval = GenPvalues[k][genid]
            GenQval = GenQvalues[k][genid]
            GPname  = 'GP_' + str(genid) + '_' + str(nodeID) + '_' + str(k)
            GQname  = 'GQ_' + str(genid) + '_' + str(nodeID) + '_' + str(k)

            GPline  = GPname + ' = ' + str(GenPval) + '\n'
            thefilevars.write(GPline)

            GQline  = GQname + ' = ' + str(GenQval) + '\n'
            thefilevars.write(GQline)

        
    log.joint(' done writing AC QCQP allvars solution to .sol file\n\n')
        
    thefilevars.close()

def writesol_k(log,all_data,k):

    casename     = all_data['casename']    
    casetype     = all_data['casetype']
    branches     = all_data['branches']
    buses        = all_data['buses']
    gens         = all_data['gens']
    IDtoCountmap = all_data['IDtoCountmap']
    objvalue     = all_data['objvalue'][k]
    vvalues      = all_data['vvalues'][k]
    thetavalues  = all_data['thetavalues'][k]
    GenPvalues   = all_data['GenPvalues'][k]
    GenQvalues   = all_data['GenQvalues'][k]
    Pfvalues     = all_data['Pfvalues'][k]
    Ptvalues     = all_data['Ptvalues'][k]
    Qfvalues     = all_data['Qfvalues'][k]
    Qtvalues     = all_data['Qtvalues'][k]
    T            = all_data['T']

    filename      = all_data['sols'] + 'ACsol_' + casename + '_' + str(k) + '_' + str(T) + '_' + casetype + '.txt'
    thefile       = open(filename,'w+')

    log.joint('\n writing solution to ' + filename + '\n')

    machinename = platform.node()
    now = time.time()
    AMPL_version = 'Version 20190223'
    solver_version = 'Artelys Knitro 14.1.0'
    opsystem = "{} {} ({})".format(platform.system(), platform.release(), platform.version())
    processor = platform.processor()
    physical_cores = psutil.cpu_count(logical=False)
    logical_processors = psutil.cpu_count(logical=True)
    cores = f"{physical_cores} physical cores, {logical_processors} logical processors"
    ram = f"{round(psutil.virtual_memory().total / (1024 ** 3))} GB RAM"
        
    thefile.write('/ACsolution : ' + casename + " " + str(k) + "/" + str(T) + " " + casetype + '\n')
    thefile.write('/Date : ' + str(time.strftime('%m-%d-%Y %H:%M:%S %Z', time.localtime(now))) + '\n')
    thefile.write('/MachineName : ' + machinename + '\n')
    thefile.write('/Processor : ' + processor + '\n')
    thefile.write('/OS : ' + opsystem + '\n')
    thefile.write('/Cores : ' + cores + '\n')
    thefile.write('/RAM : ' + ram + '\n')
    thefile.write('/AMPL : ' + AMPL_version + '\n')
    thefile.write('/Solver : ' + solver_version + '\n')
    thefile.write('objvalue ' + str(all_data['objvalue'][k]) + '\n')
    
    thefile.write('voltages and angles:\n')

    for buscount in buses.keys():
        vval       = vvalues[buscount]
        f          = buses[buscount].nodeID
        thetaval   = thetavalues[buscount]
        line = 'bus ' + str(buscount) + ' M ' + str(vval) + ' A ' + str(thetaval) + '\n'
        thefile.write(line)
            
    thefile.write('power flows and cs variables:\n')
    
    for branchid in branches.keys():

        branch     = branches[branchid]
        f          = branch.f
        t          = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        Pfval      = Pfvalues[branchid]
        Ptval      = Ptvalues[branchid]
        Qfval      = Qfvalues[branchid]
        Qtval      = Qtvalues[branchid]
        thetaval   = thetavalues[count_of_f] - thetavalues[count_of_t]
        cftval     = vvalues[count_of_f] * vvalues[count_of_t] * math.cos(thetaval)
        sftval     = vvalues[count_of_f] * vvalues[count_of_t] * math.sin(thetaval)        
        
        line = 'branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' Pft ' + str(Pfval) + ' Ptf ' + str(Ptval) + ' Qft ' + str(Qfval) + ' Qtf ' + str(Qtval) + ' cft ' + str(cftval) + ' sft ' + str(sftval) + '\n'
        thefile.write(line)

    thefile.write('generation:\n')
        
    for genid in gens.keys():
        gen     = gens[genid] 
        nodeID  = gen.nodeID
        line_gen = 'genid ' + str(genid) + ' bus ' + str(nodeID) + ' GP ' + str(GenPvalues[genid]) + ' GQ ' + str(GenQvalues[genid]) + '\n'
        thefile.write(line_gen)
        
        
    thefile.close()

    log.joint(' done writing knitro solution to .txt file\n')


def writesol_qcqp_allvars_k(log,all_data,k):

    casename      = all_data['casename']    
    casetype      = all_data['casetype']    
    branches      = all_data['branches']
    buses         = all_data['buses']
    gens          = all_data['gens']
    IDtoCountmap  = all_data['IDtoCountmap']
    objvalue      = all_data['objvalue'][k]
    vvalues       = all_data['vvalues'][k]
    thetavalues   = all_data['thetavalues'][k]
    GenPvalues    = all_data['GenPvalues'][k]
    GenQvalues    = all_data['GenQvalues'][k]
    Pfvalues      = all_data['Pfvalues'][k]
    Ptvalues      = all_data['Ptvalues'][k]
    Qfvalues      = all_data['Qfvalues'][k]
    Qtvalues      = all_data['Qtvalues'][k]
    T             = all_data['T']
    
    filenamevars  = all_data['sols'] + 'ACsol_' + casename + '_' + str(k) + '_' + str(T) + '_' + casetype + '.sol'
    thefilevars   = open(filenamevars,'w+')
    
    log.joint('\n writing solution to ' + filenamevars + '\n')


    machinename = platform.node()
    now = time.time()
    AMPL_version = 'Version 20190223'
    solver_version = 'Artelys Knitro 14.1.0'
    opsystem = "{} {} ({})".format(platform.system(), platform.release(), platform.version())
    processor = platform.processor()
    physical_cores = psutil.cpu_count(logical=False)
    logical_processors = psutil.cpu_count(logical=True)
    cores = f"{physical_cores} physical cores, {logical_processors} logical processors"
    ram = f"{round(psutil.virtual_memory().total / (1024 ** 3))} GB RAM"
    
    thefilevars.write('/ACsolution : ' + casename + " " + str(k) + "/" + str(T) + " " + casetype + '\n')    
    thefilevars.write('/Date : ' + str(time.strftime('%m-%d-%Y %H:%M:%S %Z', time.localtime(now))) + '\n')
    thefilevars.write('/MachineName : ' + machinename + '\n')
    thefilevars.write('/Processor : ' + processor + '\n')
    thefilevars.write('/OS : ' + opsystem + '\n')
    thefilevars.write('/Cores : ' + cores + '\n')
    thefilevars.write('/RAM : ' + ram + '\n')
    thefilevars.write('/AMPL : ' + AMPL_version + '\n')
    thefilevars.write('/Solver : ' + solver_version + '\n')
    thefilevars.write('/Objvalue : ' + str(all_data['objvalue'][k]) + '\n')
    
    for buscount in buses.keys():
        bus        = buses[buscount]
        vval       = vvalues[buscount]
        f          = bus.nodeID
        thetaval   = thetavalues[buscount]
        v2value    = vval**2
        evalue     = vval * math.cos(thetaval)
        fvalue     = vval * math.sin(thetaval)

        v2name     = 'c_' + str(f) + '_' + str(f)
        ename      = 'e_' + str(f)
        fname      = 'f_' + str(f)

        v2line     = v2name + ' = ' + str(v2value) + '\n'
        thefilevars.write(v2line)

        eline     = ename + ' = ' + str(evalue) + '\n'
        thefilevars.write(eline)

        fline     = fname + ' = ' + str(fvalue) + '\n'
        thefilevars.write(fline)        

        IPvalue = - bus.Pd
        IQvalue = - bus.Qd
        IPname  = 'IP_' + str(bus.nodeID)
        IQname  = 'IQ_' + str(bus.nodeID)

        for gencounter in bus.genidsbycount:
            if gens[gencounter].status:
                IPvalue += GenPvalues[gencounter]
                IQvalue += GenQvalues[gencounter]

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
        thetaftval = thetavalues[count_of_f] - thetavalues[count_of_t]
        Pfval      = Pfvalues[branchid]
        Ptval      = Ptvalues[branchid]
        Qfval      = Qfvalues[branchid]
        Qtval      = Qtvalues[branchid]
        cftval     = vvalues[count_of_f] * vvalues[count_of_t] * math.cos(thetaftval)                            
        sftval     = vvalues[count_of_f] * vvalues[count_of_t] * math.sin(thetaftval)                            

        Pfname  = 'P_' + str(branchid) + '_' + str(f) + '_' + str(t)
        Ptname  = 'P_' + str(branchid) + '_' + str(t) + '_' + str(f)
        Qfname  = 'Q_' + str(branchid) + '_' + str(f) + '_' + str(t)
        Qtname  = 'Q_' + str(branchid) + '_' + str(t) + '_' + str(f)
        cftname = 'c_' + str(branchid) + '_' + str(f) + '_' + str(t)
        sftname = 's_' + str(branchid) + '_' + str(f) + '_' + str(t)        

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

        
    for genid in gens.keys():
        gen     = gens[genid] 
        nodeID  = gen.nodeID
        GenPval = GenPvalues[genid]
        GenQval = GenQvalues[genid]
        GPname  = "GP_" + str(genid) + "_" + str(nodeID)
        GQname  = "GQ_" + str(genid) + "_" + str(nodeID)

        GPline  = GPname + ' = ' + str(GenPval) + '\n'
        thefilevars.write(GPline)

        GQline  = GQname + ' = ' + str(GenQval) + '\n'
        thefilevars.write(GQline)

        
    log.joint(' done writing AC QCQP allvars solution to .sol file\n\n')
        
    thefilevars.close()


def getloads(log,all_data):

    casename      = all_data['casename']
    T             = all_data['T']
    loadsfilename = all_data['loadsfilename']
    
    try:
        thefile = open(loadsfilename, "r")
        lines = thefile.readlines()
        lenlines = len(lines) - 1 # END                                                  
        thefile.close()
    except:
        log.stateandquit("Cannot open file " + filename)
        sys.exit("Check file with loads")

    buses        = all_data['buses']
    IDtoCountmap = all_data['IDtoCountmap']

    Pd_mtp     = {}
    for k in range(T):
        Pd_mtp[k] = {}
        
    linenum = 0

    log.joint(' reading file with multi-period demands\n')
    while linenum < lenlines:
        thisline         = lines[linenum].split()
        buscount         = int(thisline[1])
        k                = int(thisline[5])
        load             = float(thisline[7]) 
        Pd_mtp[k][buscount]  = load
        linenum         += 1

    log.joint(' done reading multi-period demands\n')        
    return  Pd_mtp


def getrampr(log,all_data):

    casename = all_data['casename']
    T        = all_data['T']
    filename = all_data['rampfilename']
    try:
        thefile = open(filename, "r")
        lines = thefile.readlines()
        lenlines = len(lines) - 1 # END                                                      
        thefile.close()
    except:
        log.stateandquit("Cannot open file " + datafilename)
        sys.exit("Check file with ramp rates")

    buses        = all_data['buses']
    IDtoCountmap = all_data['IDtoCountmap']

    rampru_mtp   = {}
    ramprd_mtp   = {}

    for k in range(T):
        rampru_mtp[k] = {}
        ramprd_mtp[k] = {}
        
    linenum = 0

    log.joint(' reading file with multi-period ramp rates\n')
    while linenum < lenlines:
        thisline = lines[linenum].split()
        gencount        = int(thisline[1])
        k               = int(thisline[5])
        rpru            = float(thisline[7])
        rprd            = float(thisline[9])
        rampru_mtp[k][gencount] = rpru
        ramprd_mtp[k][gencount] = rprd
        linenum        += 1

    log.joint(' done reading multi-period ramp rates\n')
    return  rampru_mtp, ramprd_mtp


def getsol_ac_mtp2(log,all_data): # AC one solution at a time

    T            = all_data['T']
    casename     = all_data['casename']
    casetype     = all_data['casetype']
    branches     = all_data['branches']
    buses        = all_data['buses']
    IDtoCountmap = all_data['IDtoCountmap']
    gens         = all_data['gens']

    Vinit_mtp         = {}
    thetadiffinit_mtp = {}

    
    sol_obj        = {}
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
        filename = 'ACsol_' + all_data['casename'] + '_' + str(k) + '_' + str(T) + '_' + casetype + '.txt'

        try:
            thefile   = open(filename, "r")
            lines     = thefile.readlines()
            lenlines  = len(lines)
            thefile.close()
        except:
            log.stateandquit('cannot open file ' + filename + ' k ' + str(k))
            sys.exit("failure")

        Vinit_mtp[k]          = {}
        thetadiffinit_mtp[k]  = {}
            
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
        log.joint(' reading file ' + filename + '\n')
        while linenum < lenlines:
            thisline = lines[linenum].split()
            log.joint(' the line ' + str(thisline) + '\n')
            if thisline[0] == 'objvalue':
                sol_obj[k]     = float(thisline[1])
            elif thisline[0] == 'bus':
                buscount             = int(thisline[1])
                bus                  = buses[buscount]
                sol_vm[k][bus]       = float(thisline[3])
                sol_angle[k][bus]    = float(thisline[5])
                sol_cvalues[k][bus]  = float(thisline[3])**2
                if all_data['initial_solution_0']:
                    if k == 0:
                        Vinit_mtp[k][buscount] = float(thisline[3])
                    else:
                        Vinit_mtp[k][buscount] = Vinit_mtp[buscount,0]
                else:
                    Vinit_mtp[k][buscount]    = float(thisline[3])
            elif thisline[0] == 'branch':
                branchcount             = int(thisline[1])
                branch                  = branches[branchcount]
                sol_Pfvalues[k][branch] = float(thisline[7])
                sol_Ptvalues[k][branch] = float(thisline[9])
                sol_Qfvalues[k][branch] = float(thisline[11])
                sol_Qtvalues[k][branch] = float(thisline[13])
                sol_cvalues[k][branch]  = float(thisline[15])
                sol_svalues[k][branch]  = float(thisline[17])
                if all_data['initial_solution_0']:
                    if k == 0:
                        thetadiffinit_mtp[k][branchcount] = math.atan2(float(thisline[17]),float(thisline[15]))
                    else:
                        thetadiffinit_mtp[k][branchcount] = thetadiffinit_mtp[0][branchcoun]
                else:
                    thetadiffinit_mtp[k][branchcount] = math.atan2(float(thisline[17]),float(thisline[15]))
            elif thisline[0] == 'genid':
                genid      = int(thisline[1])
                gen_nodeID = int(thisline[3])
                sol_GenPvalues[k][genid] = float(thisline[5])
                sol_GenQvalues[k][genid] = float(thisline[7])      

            linenum += 1

        log.joint(' done loading solution to ' + str(k) + 'th time-period\n')

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

    log.joint(' knitroampl multi-time period solution loaded\n')

    return Vinit_mtp, thetadiffinit_mtp


def initial_solution_mtp2(log,all_data):

    T = all_data['T']
    
    if all_data['initial_solution'] and all_data['ampl_sol']:
        log.joint(' loading initial solution...\n')
        Vinit_mtp, thetadiffinit_mtp = getsol_ac_mtp2(log,all_data) # Using AC solutions obtained a period at a time
    else:
        log.joint(' flat-start (default)\n')
        Vinit_mtp         = {}
        thetadiffinit_mtp = {}
        for k in range(T):
            Vinit_mtp[k]         = {}
            thetadiffinit_mtp[k] = {}
            
            for buscount in all_data['buses'].keys():
                Vinit_mtp[k][buscount] = 1  

            for branchcount in all_data['branches'].keys():
                thetadiffinit_mtp[k][branchcount] = 0

    return Vinit_mtp, thetadiffinit_mtp


def datastruct_acsolutions(log,all_data):

    T = all_data['T']
    
    all_data['objvalue']      = {}
    all_data['vvalues']       = {}
    all_data['thetavalues']   = {}
    all_data['GenPvalues']    = {}
    all_data['GenQvalues']    = {}
    all_data['Pfvalues']      = {}
    all_data['Ptvalues']      = {}
    all_data['Qfvalues']      = {}
    all_data['Qtvalues']      = {}
    all_data['sumPd']         = {}
    all_data['sumQd']         = {}
    all_data['sumPLoss']      = {}
    all_data['sumQLoss']      = {}
    all_data['sumGenP']       = {}
    all_data['sumGenQ']       = {}    
    all_data['solver_status'] = {}
    
    for k in range(T):
        all_data['objvalue'][k]    = {}
        all_data['vvalues'][k]     = {}
        all_data['thetavalues'][k] = {}
        all_data['GenPvalues'][k]  = {}
        all_data['GenQvalues'][k]  = {}
        all_data['Pfvalues'][k]    = {}
        all_data['Ptvalues'][k]    = {}
        all_data['Qfvalues'][k]    = {}
        all_data['Qtvalues'][k]    = {}


def datastruct_Pmax_min(log,all_data):

    T = all_data['T']
    
    Pmax_mtp = {}
    Pmin_mtp = {}
    
    for k in range(T):
        Pmax_mtp[k] = {}
        Pmin_mtp[k] = {}

    return Pmax_mtp, Pmin_mtp
