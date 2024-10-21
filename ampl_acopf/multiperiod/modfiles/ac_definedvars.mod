###############################################################################
##                                                                           ##      
## This AMPL modfile of an ACOPF polar formulation was written and           ##
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
param theta_max {i in buses};
param theta_min {i in buses};
param maxangle {i in branches};
param minangle {i in branches};
param Vinit {i in buses}; #initial point
param thetadiffinit {i in branches}; #initial point			    

#VARIABLES
					   
var v {i in buses} >= Vmin[i], <= Vmax[i], := Vinit[i]; 
var theta {i in buses} >= theta_min[i], <= theta_max[i];
var thetadiff {i in branches} >= minangle[i], <= maxangle[i], := thetadiffinit[i];
#var Pf {i in branches} >= - U[i], <= U[i]; 
#var Pt {i in branches} >= - U[i], <= U[i]; 
#var Qf {i in branches} >= - U[i], <= U[i]; 
#var Qt {i in branches} >= - U[i], <= U[i]; 
var Pg {i in gens} >= Pmin[i], <= Pmax[i]; 
var Qg {i in gens} >= Qmin[i], <= Qmax[i]; 


#DEFINED VARS

#defP

var Pf{i in branches} =
                   Gff[i] * v[bus_f[i]] * v[bus_f[i]] + Gft[i] * v[bus_f[i]] * v[bus_t[i]] * cos(thetadiff[i]) + Bft[i] * v[bus_f[i]] * v[bus_t[i]] * sin(thetadiff[i]);

var Pt{i in branches} =
                   Gtt[i] * v[bus_t[i]] * v[bus_t[i]] + Gtf[i] * v[bus_f[i]] * v[bus_t[i]] * cos(thetadiff[i]) - Btf[i] * v[bus_f[i]] * v[bus_t[i]] * sin(thetadiff[i]);

#defQ

var Qf{i in branches} =
                   - Bff[i] * v[bus_f[i]] * v[bus_f[i]] - Bft[i] * v[bus_f[i]] * v[bus_t[i]] * cos(thetadiff[i]) + Gft[i] * v[bus_f[i]] * v[bus_t[i]] * sin(thetadiff[i]);

var Qt{i in branches} =
                   - Btt[i] * v[bus_t[i]] * v[bus_t[i]] - Btf[i] * v[bus_f[i]] * v[bus_t[i]] * cos(thetadiff[i]) - Gtf[i] * v[bus_f[i]] * v[bus_t[i]] * sin(thetadiff[i]);

#OBJECTIVE                                                  

minimize total_cost: sum {i in gens} (fixedcost[i] + lincost[i] * Pg[i] + quadcost[i] * Pg[i] ^2);
           

#power balance

subject to Pbalance {i in buses}:
       (if i in bus_Gs then (Gs[i] * v[i] * v[i] ) else 0  ) + (sum {j in branches_f[i]} Pf[j] ) + (sum {j in branches_t[i]} Pt[j] ) = (if card(bus_gens[i]) > 0 then (sum {k in bus_gens[i]} Pg[k] ) else 0) - Pd[i];

subject to Qbalance {i in buses}:
       (if i in bus_Bs then (- Bs[i] * v[i] * v[i] ) else 0  ) + (sum {j in branches_f[i]} Qf[j] ) + (sum {j in branches_t[i]} Qt[j] ) = (if card(bus_gens[i]) > 0 then (sum {k in bus_gens[i]} Qg[k] ) else 0) - Qd[i];


#limits

subject to limits_f {i in branches}: Pf[i] ^2 + Qf[i] ^2 <= U[i] ^2;
subject to limits_t {i in branches}: Pt[i] ^2 + Qt[i] ^2 <= U[i] ^2;

#angle diff

subject to thetadiff_def {i in branches}: thetadiff[i] = theta[bus_f[i]] - theta[bus_t[i]];



             
