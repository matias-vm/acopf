###############################################################################
##                                                                           ##
## This code was written and is being maintained by Matias Villagra,         ##
## PhD Student in Operations Research @ Columbia, supervised by              ##
## Daniel Bienstock.                                                         ##
##                                                                           ##
## Please report any bugs or issues (for sure there will be) to              ##
##                        mjv2153@columbia.edu                               ##
##                                                                           ##
## Oct 2023                                                                  ##
###############################################################################

from amplpy import AMPL
import numpy as np
from myutils import breakexit
from log import danoLogger
import time
import math
import os
import platform
import psutil

def goac(log,all_data):

    log.joint(' creating ampl object ...\n')
    ampl                    = AMPL()
    all_data['ampl_object'] = ampl

    IDtoCountmap = all_data['IDtoCountmap']    
    modfile      = all_data['modfile']
    solver       = all_data['solver']

    log.joint(" reading modfile ...\n")
    
    t0 = time.time()
    ampl.read('../modfiles/'+modfile)
    t1 = time.time()
    log.joint(" modfile read in time " + str(t1-t0))

    ampl.eval("option display_precision 0;")
    ampl.eval("option expand_precision 0;")
    ampl.setOption('solver',solver)
    
    ampl.setOption('presolve',0)
    log.joint(' AMPL presolve off\n')
        
    log.joint(" solver set to " + solver + "\n")
    
    ampl.eval("option knitro_options 'feastol_abs=1e-6 opttol_abs=1e-6 blasoptionlib=1 numthreads=20 linsolver=5 maxtime_real=" + str(all_data['max_time']) + "';")

                
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
    thetainit    = {}


    for bus in all_data['buses'].values():
        buscount = bus.count

        if ( len(bus.frombranchids) == 0 ) and ( len(bus.tobranchids) == 0 ): 
            log.joint(' isolated bus ' + str(bus.nodeID) + ' busid ' + str(buscount) + '\n')
        
        buses[buscount] = bus.nodeID 
        Pd[buscount] = bus.Pd
        Qd[buscount] = bus.Qd

        Vmin[buscount] = bus.Vmin
        Vmax[buscount] = bus.Vmax
        
        Vinit[buscount]      = 1
        thetainit[buscount]  = 0

        branches_f[buscount] = []
        branches_t[buscount] = []

        bus_gens[buscount] = bus.genidsbycount 

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
    ampl.get_parameter('Pd').set_values(Pd)
    ampl.get_parameter('Qd').set_values(Qd)
    ampl.get_parameter('Vmax').set_values(Vmax)
    ampl.get_parameter('Vmin').set_values(Vmin)
    ampl.get_parameter('theta_min').setValues(theta_min)
    ampl.get_parameter('theta_max').setValues(theta_max)

    ampl.get_parameter('Vinit').set_values(Vinit)
    #ampl.get_parameter('thetainit').set_values(thetainit)

    branches_f_set = ampl.getSet('branches_f')
    branches_t_set = ampl.getSet('branches_t')
    bus_gens_set   = ampl.getSet('bus_gens')
    
    for bus in all_data['buses'].values():
        branches_f_set[bus.count].setValues(branches_f[bus.count])
        branches_t_set[bus.count].setValues(branches_t[bus.count])
        bus_gens_set[bus.count].setValues(bus_gens[bus.count])
    
    # Setting up branches
    branches  = {}
    U         = {}
    Gtt       = {}
    Btt       = {}
    Gff       = {}
    Bff       = {}
    Gtf       = {}
    Btf       = {}
    Gft       = {}
    Bft       = {}
    bus_f     = {}
    bus_t     = {}
    maxangle  = {}
    minangle  = {}
    thetadiffinit = {}
    
    for branch in all_data['branches'].values():
        branchcount            = branch.count

        if branch.status != 1: #the reader does not consider branches whose status is 0, hence this should never be true
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
        
        thetadiffinit[branchcount]  = 0   
            
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
    ampl.get_parameter('thetadiffinit').setValues(thetadiffinit)

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
        Pmax[gencount]      = gen.Pmax * gen.status
        Pmin[gencount]      = gen.Pmin * gen.status
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
    ampl.get_parameter('Pmax').setValues(Pmax)
    ampl.get_parameter('Pmin').setValues(Pmin)
    ampl.get_parameter('Qmax').setValues(Qmax)
    ampl.get_parameter('Qmin').setValues(Qmin)
    ampl.get_parameter('fixedcost').setValues(fixedcost)
    ampl.get_parameter('lincost').setValues(lincost)
    ampl.get_parameter('quadcost').setValues(quadcost)
    
    log.joint(" sets and parameters loaded\n")

    log.joint(" saving processed data to all_data\n")

    expand = all_data['expand']

    # Expands the model and write a human-readable file (similar in spirit to an .lp file)
    if expand:
        time1 = time.time()
        filename = 'basemodel.out'
        doNLP = 1
        if doNLP:
            filename = 'NLP.out'
        log.joint('Now expanding to %s.\n'%(filename))
        amplstate = 'expand; display {j in 1.._nvars} (_varname[j],_var[j].lb,_var[j].ub);'
        modelout = ampl.getOutput(amplstate)
        outfile = open(filename,"w")
        outfile.write("model = " + str(modelout) + "\n")
        outfile.close()
        
    # Solve
    log.joint(" solving model ...\n")
    t0 = time.time()
    ampl.solve()
    t1 = time.time()


    log.joint("\n ===============================================================\n")
    log.joint(" ===============================================================\n")
    
    # Get solution
    total_costvar           = ampl.get_objective("total_cost")
    vvar                    = ampl.get_variable("v")
    thetavar                = ampl.get_variable("theta")
    Pfvar                   = ampl.get_variable("Pf")
    Ptvar                   = ampl.get_variable("Pt")
    Qfvar                   = ampl.get_variable("Qf")
    Qtvar                   = ampl.get_variable("Qt")
    GenPvar                 = ampl.get_variable("Pg")
    GenQvar                 = ampl.get_variable("Qg")
    all_data['objvalue']    = total_costvar.get().value()
    all_data['vvalues']     = vvar.get_values().to_dict()
    all_data['thetavalues'] = thetavar.get_values().to_dict()
    all_data['GenPvalues']  = GenPvar.get_values().to_dict()
    all_data['GenQvalues']  = GenQvar.get_values().to_dict()
    all_data['Pfvalues']    = Pfvar.get_values().to_dict()
    all_data['Ptvalues']    = Ptvar.get_values().to_dict()
    all_data['Qfvalues']    = Qfvar.get_values().to_dict()
    all_data['Qtvalues']    = Qtvar.get_values().to_dict()

    
    PLoss   = {}
    QLoss   = {}
    for branch in all_data['Pfvalues'].keys():
        PLoss[branch] = all_data['Pfvalues'][branch] + all_data['Ptvalues'][branch] 
        QLoss[branch] = all_data['Qfvalues'][branch] + all_data['Qtvalues'][branch] 

    timesofar = time.time() - all_data['T0']        
        
    log.joint(" case " + all_data['casefilename'] + "\n")
    log.joint(" modfile " + all_data['modfile'] + "\n")
    log.joint(" solver " + all_data['solver'] + "\n")
    log.joint(" objective " + str(all_data['objvalue']) + '\n')
    log.joint(" active power generation " + str(sum(all_data['GenPvalues'])) + '\n')
    log.joint(" active power demand " + str(sum(Pd.values())) + '\n')
    log.joint(" active power loss " + str(sum(PLoss.values())) + '\n')
    log.joint(" reactive power generation " + str(sum(all_data['GenQvalues'])) + '\n')
    log.joint(" reactive power demand " + str(sum(Qd.values())) + '\n')
    log.joint(" reactive power loss " + str(sum(QLoss.values())) + '\n')
    log.joint(" solver runtime " + str(t1-t0) + '\n')
    log.joint(" time so far " + str(timesofar) + '\n')
    
    log.joint(" writing casename, modfile, obj and runtime to summary_ac.log\n")

    summary_ac = open("summary_ac.log","a+")

    summary_ac.write(' case ' + all_data['casename'] + ' modfile ' + all_data['modfile']
                     + ' solver ' + all_data['solver'] + ' obj ' + str(all_data['objvalue'])
                     + ' runtime ' + str(timesofar) + '\n')

    summary_ac.close()

    log.joint(" ===============================================================\n")
    log.joint(" ===============================================================\n")
    
    if all_data['writesol']:
        writesol(log,all_data)
        writesol_qcqp_allvars(log,all_data)

            
def writesol(log,all_data):

    ampl         = all_data['ampl_object']
    branches     = all_data['branches']
    buses        = all_data['buses']
    gens         = all_data['gens']
    IDtoCountmap = all_data['IDtoCountmap']
    tolerance    = 1e-05

    objvalue     = all_data['objvalue']
    vvalues      = all_data['vvalues']
    thetavalues  = all_data['thetavalues']
    GenPvalues   = all_data['GenPvalues']
    GenQvalues   = all_data['GenQvalues']
    Pfvalues     = all_data['Pfvalues']
    Ptvalues     = all_data['Ptvalues']
    Qfvalues     = all_data['Qfvalues']
    Qtvalues     = all_data['Qtvalues']
    
    filename      = 'ACsol_' + all_data['casename'] + '.txt'
    thefile       = open(filename,'w+')

    log.joint(' writing solution to ' + filename + '\n')

    # Get machine name and current time
    machinename = platform.node()
    now         = time.time()

    # Get AMPL version and solver
    version_str    = ampl.get_option('version')
    AMPL_version   = f"version date {version_str.split()[2]}"
    solver_version = all_data['solver']

    # System information
    opsystem           = f"{platform.system()} {platform.release()} ({platform.version()})"
    processor          = platform.processor()
    physical_cores     = psutil.cpu_count(logical=False)
    logical_processors = psutil.cpu_count(logical=True)
    cores              = f"{physical_cores} physical cores, {logical_processors} logical processors"
    ram                = f"{round(psutil.virtual_memory().total / (1024 ** 3))} GB RAM"
    
    thefile.write('/ACsolution : ' + all_data['casename'] + '\n')
    thefile.write('/Date : ' + str(time.strftime('%m-%d-%Y %H:%M:%S %Z', time.localtime(now))) + '\n')
    thefile.write('/MachineName : ' + machinename + '\n')
    thefile.write('/Processor : ' + processor + '\n')
    thefile.write('/OS : ' + opsystem + '\n')
    thefile.write('/Cores : ' + cores + '\n')
    thefile.write('/RAM : ' + ram + '\n')
    thefile.write('/AMPL : ' + AMPL_version + '\n')
    thefile.write('/Solver : ' + solver_version + '\n')
    thefile.write('objvalue ' + str(all_data['objvalue']) + '\n')
    
    thefile.write('voltages and angles:\n')

    for buscount in buses.keys():
        vval       = vvalues[buscount]
        f          = buses[buscount].nodeID
        thetaval   = thetavalues[buscount] * (180 / math.pi)  # angles in degrees
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

    log.joint(' done writing knitro solution to .txt file\n\n')

def writesol_qcqp_allvars(log,all_data):

    ampl          = all_data['ampl_object']
    branches      = all_data['branches']
    buses         = all_data['buses']
    gens          = all_data['gens']
    IDtoCountmap  = all_data['IDtoCountmap']
    tolerance     = 1e-05

    objvalue     = all_data['objvalue']
    vvalues      = all_data['vvalues']
    thetavalues  = all_data['thetavalues']
    GenPvalues   = all_data['GenPvalues']
    GenQvalues   = all_data['GenQvalues']
    Pfvalues     = all_data['Pfvalues']
    Ptvalues     = all_data['Ptvalues']
    Qfvalues     = all_data['Qfvalues']
    Qtvalues     = all_data['Qtvalues']


    filenamevars  = 'ACsol_' + all_data['casename'] + '.sol'    
    thefilevars   = open(filenamevars,'w+')
    
    log.joint(' writing solution to ' + filenamevars + '\n')

    # Get machine name and current time
    machinename = platform.node()
    now         = time.time()

    # Get AMPL version and solver
    version_str    = ampl.get_option('version')
    AMPL_version   = f"version date {version_str.split()[2]}"
    solver_version = all_data['solver']

    # System information
    opsystem           = f"{platform.system()} {platform.release()} ({platform.version()})"
    processor          = platform.processor()
    physical_cores     = psutil.cpu_count(logical=False)
    logical_processors = psutil.cpu_count(logical=True)
    cores              = f"{physical_cores} physical cores, {logical_processors} logical processors"
    ram                = f"{round(psutil.virtual_memory().total / (1024 ** 3))} GB RAM"
    
    thefilevars.write('/ACsolution : ' + all_data['casename'] + '\n')
    thefilevars.write('/Date : ' + str(time.strftime('%m-%d-%Y %H:%M:%S %Z', time.localtime(now))) + '\n')
    thefilevars.write('/MachineName : ' + machinename + '\n')
    thefilevars.write('/Processor : ' + processor + '\n')
    thefilevars.write('/OS : ' + opsystem + '\n')
    thefilevars.write('/Cores : ' + cores + '\n')
    thefilevars.write('/RAM : ' + ram + '\n')
    thefilevars.write('/AMPL : ' + AMPL_version + '\n')
    thefilevars.write('/Solver : ' + solver_version + '\n')
    thefilevars.write('/Objvalue : ' + str(all_data['objvalue']) + '\n')
    
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
