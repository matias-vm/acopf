###############################################################################
##                                                                           ##      
## This AMPL modfile of the Jabr relaxation for ACOPF was written and        ##
## is being maintained by Matias Villagra,                                   ##
## PhD Student in Operations Research @ Columbia, supervised by              ## 
## Daniel Bienstock.                                                         ##      
##                                                                           ##    
## Please report any bugs or issues to mjv2153@columbia.edu                  ##     
##                                                                           ##      
## Oct 2023                                                                  ## 
###############################################################################


#SETS                                
set buses;
set gens;
set branches;
set branches_f {i in buses};
set branches_t {i in buses};
set bus_gens {i in buses};
set bus_Gs;
set bus_Bs;

#PARAMETERS
param fixedcost {i in gens};
param lincost {i in gens}; 
param quadcost {i in gens}; 

param Gtt {i in branches}; 
param Btt {i in branches};
param Gff {i in branches}; 
param Bff {i in branches}; 
param Gtf {i in branches}; 
param Btf {i in branches}; 
param Gft {i in branches}; 
param Bft {i in branches}; 
param Rft {i in branches};

param bus_f {i in branches};
param bus_t {i in branches};    
param Gs {i in bus_Gs};
param Bs {i in bus_Bs};
                                                                              
param Vmax {i in buses} >= 0; #max voltage (non squared)
param Vmin {i in buses} >= 0; #min voltage (non squared)
param Pmax {i in gens}; #max active power gen
param Pmin {i in gens}; #min active power gen
param Qmax {i in gens}; #max reactive power gen
param Qmin {i in gens}; #min reactive power gen                            
param Pd {i in buses}; #active power demand
param Qd {i in buses}; #reactive power demand
param U {i in branches} >= 0; #branch power flow limit
param Vinit {i in buses};
param Cinit {i in branches};
param Sinit {i in branches};
param c_ubound {i in branches};
param c_lbound {i in branches};
param s_ubound {i in branches};
param s_lbound {i in branches};

#VARIABLES                                
var c {i in branches} >= c_lbound[i], <= c_ubound[i];#, := Cinit[i];
var s {i in branches} >= s_lbound[i], <= s_ubound[i];#, := Sinit[i]; 
var v {i in buses} >= Vmin[i] ^2, <= Vmax[i] ^2, := Vinit[i];
var Pf {i in branches} >= - U[i], <= U[i]; #limit branches
var Pt {i in branches} >= - U[i], <= U[i]; #limit branches
var Qf {i in branches} >= - U[i], <= U[i]; #limit branches
var Qt {i in branches} >= - U[i], <= U[i]; #limit branches                   
var Pg {i in gens} >= Pmin[i], <= Pmax[i]; #active power generator
var Qg {i in gens} >= Qmin[i], <= Qmax[i]; #reactive power generator          

#OBJECTIVE                                                  

#minimize total_cost: sum {i in gens} Pg[i];
minimize total_cost: sum {i in gens} (fixedcost[i] + lincost[i] * Pg[i] + quadcost[i] * Pg[i] ^2);
           
#CONSTRAINTS

#def PandQ

subject to Pf_def {i in branches}:
                   Pf[i] = Gff[i] * v[bus_f[i]] + Gft[i] * c[i] + Bft[i] * s[i];

subject to Pt_def {i in branches}: 
                   Pt[i] = Gtt[i] * v[bus_t[i]] + Gtf[i] * c[i] - Btf[i] * s[i];

subject to Qf_def {i in branches}:
                   Qf[i] = - Bff[i] * v[bus_f[i]] - Bft[i] * c[i] + Gft[i] * s[i];

subject to Qt_def {i in branches}:
                   Qt[i] = - Btt[i] * v[bus_t[i]] - Btf[i] * c[i] - Gtf[i] * s[i];

#power balance

subject to Pbalance {i in buses}:
       (if i in bus_Gs then (Gs[i] * v[i] ) else 0  ) + (sum {j in branches_f[i]} Pf[j] ) + (sum {j in branches_t[i]} Pt[j] ) = (if card(bus_gens[i]) > 0 then (sum {k in bus_gens[i]} Pg[k] ) else 0) - Pd[i];

subject to Qbalance {i in buses}:
       (if i in bus_Bs then (- Bs[i] * v[i] ) else 0  ) + (sum {j in branches_f[i]} Qf[j] ) + (sum {j in branches_t[i]} Qt[j] ) = (if card(bus_gens[i]) > 0 then (sum {k in bus_gens[i]} Qg[k] ) else 0) - Qd[i];

#jabrs
subject to jabr {i in branches}: c[i] ^2 + s[i] ^2 <= v[bus_f[i]] * v[bus_t[i]];


#limits

subject to limits_f {i in branches}: Pf[i] ^2 + Qf[i] ^2 <= U[i] ^2;

subject to limits_t {i in branches}: Pt[i] ^2 + Qt[i] ^2 <= U[i] ^2;
             
