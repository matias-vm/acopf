###############################################################################
##                                                                           ##      
## This AMPL modfile of a conic relaxation (i2 SOCP) for ACOPF was           ##
## written and is being maintained by Matias Villagra,                       ##
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
set quadgens;
set lingens;
     
#PARAMETERS
param fixedcost {i in gens};
param lincost {i in gens}; 
param quadcost {i in gens};
param alfa {i in quadgens};
param beta {i in quadgens};
param gamma {i in quadgens};
param sigma;

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
param i2max {i in branches} >= 0;
param g {i in branches};
param b {i in branches};
param bshunt {i in branches};
param ratio {i in branches};
param phase_angle {i in branches};
param CSmax {i in branches} >= 0;
param Cinit {i in branches};
param Sinit {i in branches};
param c_ubound {i in branches};
param c_lbound {i in branches};
param s_ubound {i in branches};
param s_lbound {i in branches};

#VARIABLES                                
#var c {i in branches} >= - CSmax[i], <= CSmax[i];#, := Cinit[i];
#var s {i in branches} >= - CSmax[i], <= CSmax[i];#, := Sinit[i];
var c {i in branches} >= c_lbound[i], <= c_ubound[i];#, := Cinit[i];
var s {i in branches} >= s_lbound[i], <= s_ubound[i];#, := Sinit[i]; 
var v {i in buses} >= Vmin[i] ^2, <= Vmax[i] ^2, := Vinit[i];
var Pf {i in branches} >= - U[i], <= U[i]; #limit branches
var Pt {i in branches} >= - U[i], <= U[i]; #limit branches
var Qf {i in branches} >= - U[i], <= U[i]; #limit branches
var Qt {i in branches} >= - U[i], <= U[i]; #limit branches                 
var Pg {i in gens} >= Pmin[i], <= Pmax[i]; #active power generator
var Qg {i in gens} >= Qmin[i], <= Qmax[i]; #reactive power generator
var i2f {i in branches} >= 0, <= i2max[i]; 

var t_quadcost >=0;
var t_lincost >= 0;
var scaledPg {i in quadgens} >= 0;


#OBJECTIVE                                                  

minimize total_cost: t_quadcost + t_lincost + sigma;
           
#CONSTRAINTS

#def obj

subject to quadobj: sum {i in quadgens} scaledPg[i] ^2 <= t_quadcost;

subject to linobj: sum {i in lingens} (fixedcost[i] + lincost[i] * Pg[i]) <= t_lincost;


#def PandQ

subject to Pf_def {i in branches}:
                   Pf[i] = Gff[i] * v[bus_f[i]] + Gft[i] * c[i] + Bft[i] * s[i];

subject to Pt_def {i in branches}: 
                   Pt[i] = Gtt[i] * v[bus_t[i]] + Gtf[i] * c[i] - Btf[i] * s[i];

subject to Qf_def {i in branches}:
                   Qf[i] = - Bff[i] * v[bus_f[i]] - Bft[i] * c[i] + Gft[i] * s[i];

subject to Qt_def {i in branches}:
                   Qt[i] = - Btt[i] * v[bus_t[i]] - Btf[i] * c[i] - Gtf[i] * s[i];


#def i2 var
subject to i2f_def {i in branches}:
        i2f[i] = ( ( g[i] * g[i] + b[i] * b[i] ) / ( ratio[i] * ratio[i] ) ) * ( ( v[bus_f[i]] / (ratio[i] * ratio[i]) ) + v[bus_t[i]] - ( 2 / ratio[i] ) * ( c[i] * cos(phase_angle[i]) + s[i] * sin(phase_angle[i]) ) ) + ( b[i] * bshunt[i] / ( ratio[i] * ratio[i] * ratio[i] ) ) * ( ( v[bus_f[i]] / ratio[i] ) - ( c[i] * cos(phase_angle[i]) + s[i] * sin(phase_angle[i]) ) ) + ( g[i] * bshunt[i] / ( ratio[i] * ratio[i] * ratio[i] ) ) * ( s[i] * cos(phase_angle[i]) - c[i] * sin(phase_angle[i]) ) + (bshunt[i] * bshunt[i] * v[bus_f[i]] / ( 4 * ratio[i] * ratio[i] * ratio[i] * ratio[i] ) );


#power balance

subject to Pbalance {i in buses}:
       (if i in bus_Gs then (Gs[i] * v[i] ) else 0  ) + (sum {j in branches_f[i]} Pf[j] ) + (sum {j in branches_t[i]} Pt[j] ) = (if card(bus_gens[i]) > 0 then (sum {k in bus_gens[i]} Pg[k] ) else 0) - Pd[i];

subject to Qbalance {i in buses}:
       (if i in bus_Bs then (- Bs[i] * v[i] ) else 0  ) + (sum {j in branches_f[i]} Qf[j] ) + (sum {j in branches_t[i]} Qt[j] ) = (if card(bus_gens[i]) > 0 then (sum {k in bus_gens[i]} Qg[k] ) else 0) - Qd[i];


#i2s
      
subject to i2 {i in branches}: Pf[i] ^2 + Qf[i] ^2 <= v[bus_f[i]] * i2f[i];


#limits

subject to limits_f {i in branches}: Pf[i] ^2 + Qf[i] ^2 <= U[i] ^2;

subject to limits_t {i in branches}: Pt[i] ^2 + Qt[i] ^2 <= U[i] ^2;

#def scaledPg

subject to scaledPg_def {i in quadgens}: scaledPg[i] = alfa[i] * Pg[i] - beta[i];




