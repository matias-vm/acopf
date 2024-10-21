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

param T; # number of time periods
param fixedcost {i in gens};
param lincost {i in gens}; 
param quadcost {i in gens}; 

#branch parameters
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
param Pd {i in buses, h in 0..T-1}; #active power demand
param Qd {i in buses}; #reactive power demand
param U {i in branches} >= 0; #branch power flow limit
param Vinit {i in buses, h in 0..T-1};
param Cinit {i in branches, h in 0..T-1};
param Sinit {i in branches, h in 0..T-1};
param c_ubound {i in branches};
param c_lbound {i in branches};
param s_ubound {i in branches};
param s_lbound {i in branches};

param rampru {i in gens, h in 0..T-1};
param ramprd {i in gens, h in 0..T-1};

#VARIABLES                                
var c {i in branches, h in 0..T-1} >= c_lbound[i], <= c_ubound[i], := Cinit[i,h];
var s {i in branches, h in 0..T-1} >= s_lbound[i], <= s_ubound[i], := Sinit[i,h]; 
var v {i in buses, h in 0..T-1} >= Vmin[i] ^2, <= Vmax[i] ^2, := Vinit[i,h];
var Pf {i in branches, h in 0..T-1} >= - U[i], <= U[i]; #limit branches
var Pt {i in branches, h in 0..T-1} >= - U[i], <= U[i]; #limit branches
var Qf {i in branches, h in 0..T-1} >= - U[i], <= U[i]; #limit branches
var Qt {i in branches, h in 0..T-1} >= - U[i], <= U[i]; #limit branches             
var Pg {i in gens, h in 0..T-1} >= Pmin[i], <= Pmax[i]; #active power generator
var Qg {i in gens, h in 0..T-1} >= Qmin[i], <= Qmax[i]; #reactive power generator     

#OBJECTIVE                                                  

minimize total_cost: sum {i in gens, h in 0..T-1} (fixedcost[i] + lincost[i] * Pg[i,h] + quadcost[i] * Pg[i,h] ^2);
           
#CONSTRAINTS

#def PandQ

subject to Pf_def {i in branches, h in 0..T-1}:
                   Pf[i,h] = Gff[i] * v[bus_f[i],h] + Gft[i] * c[i,h] + Bft[i] * s[i,h];

subject to Pt_def {i in branches, h in 0..T-1}: 
                   Pt[i,h] = Gtt[i] * v[bus_t[i],h] + Gtf[i] * c[i,h] - Btf[i] * s[i,h];

subject to Qf_def {i in branches, h in 1..T-1}:
                   Qf[i,h] = - Bff[i] * v[bus_f[i],h] - Bft[i] * c[i,h] + Gft[i] * s[i,h];

subject to Qt_def {i in branches, h in 0..T-1}:
                   Qt[i,h] = - Btt[i] * v[bus_t[i],h] - Btf[i] * c[i,h] - Gtf[i] * s[i,h];

#power balance

subject to Pbalance {i in buses, h in 0..T-1}:
       (if i in bus_Gs then (Gs[i] * v[i,h] ) else 0  ) + (sum {j in branches_f[i]} Pf[j,h] ) + (sum {j in branches_t[i]} Pt[j,h] ) = (if card(bus_gens[i]) > 0 then (sum {k in bus_gens[i]} Pg[k,h] ) else 0) - Pd[i,h];

subject to Qbalance {i in buses, h in 0..T-1}:
       (if i in bus_Bs then (- Bs[i] * v[i,h] ) else 0  ) + (sum {j in branches_f[i]} Qf[j,h] ) + (sum {j in branches_t[i]} Qt[j,h] ) = (if card(bus_gens[i]) > 0 then (sum {k in bus_gens[i]} Qg[k,h] ) else 0) - Qd[i];

#jabrs
subject to jabr {i in branches, h in 0..T-1}: c[i,h] ^2 + s[i,h] ^2 <= v[bus_f[i],h] * v[bus_t[i],h];


#limits

subject to limits_f {i in branches, h in 0..T-1}: Pf[i,h] ^2 + Qf[i,h] ^2 <= U[i] ^2;

subject to limits_t {i in branches, h in 0..T-1}: Pt[i,h] ^2 + Qt[i,h] ^2 <= U[i] ^2;

#ramping constraints

subject to ramp_up {i in gens, h in 0..T-2}: Pg[i,h+1] - (1 + rampru[i,h+1]) * Pg[i,h] <= 0;

subject to ramp_down {i in gens, h in 0..T-2}: (1 - ramprd[i,h+1]) * Pg[i,h] - Pg[i,h+1] <= 0;
