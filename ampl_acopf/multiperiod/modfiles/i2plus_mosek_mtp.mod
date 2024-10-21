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
set quadgens; #set of gens with quad cost
set lingens; #set of gens with linear cost
set goodi2; #set of goodi2 branches, see Section 4.5
set badi2; #set of badi2 branches, see Section 4.5
     
#PARAMETERS
param T; #number of time periods
param fixedcost {i in gens}; #constant term objective
param lincost {i in gens};  #coeff linear terms objective
param quadcost {i in gens}; #coeff quad terms objective
param alfa {i in quadgens}; 
param beta {i in quadgens};
param gamma {i in quadgens};
param sigma;

#bad i2 coeffs, see Section 4.5
param beta_badi2 {i in badi2}; 
param gamma_badi2 {i in badi2};
param zeta_badi2 {i in badi2};
param rhs_badi2 {i in badi2};

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
param Vinit {i in buses, h in 0..T-1}; #initial point voltages (default flat-start)
param Cinit {i in branches, h in 0..T-1}; #initial point c var (default flat-start)
param Sinit {i in branches, h in 0..T-1}; #initial point s var (default flat-start)

#branch parameters of goodi2, badi2 are not added, see Section 4.5
param i2max {i in goodi2} >= 0;
param g {i in goodi2}; 
param b {i in goodi2};
param bshunt {i in goodi2};
param ratio {i in goodi2};
param phase_angle {i in goodi2};

param c_ubound {i in branches}; #c var upper bound
param c_lbound {i in branches}; #c var lower bound
param s_ubound {i in branches}; #s var upper bound
param s_lbound {i in branches}; #s var lower bound
param rampru {i in gens, h in 0..T-1}; #ramp up rates
param ramprd {i in gens, h in 0..T-1}; #ramp down rates

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
var i2f {i in goodi2, h in 0..T-1} >= 0, <= i2max[i]; 
var t_quadcost >=0;
var t_lincost >= 0;
var scaledPg {i in quadgens, h in 0..T-1} >= 0; #Pgs are rescaled/change of variables (to make modfile compatible with the 3 solvers)
var absPg {i in gens, h in 0..T-1} >= 0; #absolute value of Pg

#OBJECTIVE                                                  

minimize total_cost: t_quadcost + t_lincost + T * sigma; #eq (1a)
           
#CONSTRAINTS

#def obj

subject to quadobj: sum {i in quadgens, h in 0..T-1} scaledPg[i,h] ^2 <= t_quadcost;

subject to linobj: sum {i in lingens, h in 0..T-1} (fixedcost[i] + lincost[i] * Pg[i,h]) <= t_lincost;


#def PandQ, eqs (6a)-(6d)

subject to Pf_def {i in branches, h in 0..T-1}:
                   Pf[i,h] = Gff[i] * v[bus_f[i],h] + Gft[i] * c[i,h] + Bft[i] * s[i,h];

subject to Pt_def {i in branches, h in 0..T-1}: 
                   Pt[i,h] = Gtt[i] * v[bus_t[i],h] + Gtf[i] * c[i,h] - Btf[i] * s[i,h];

subject to Qf_def {i in branches, h in 1..T-1}:
                   Qf[i,h] = - Bff[i] * v[bus_f[i],h] - Bft[i] * c[i,h] + Gft[i] * s[i,h];

subject to Qt_def {i in branches, h in 0..T-1}:
                   Qt[i,h] = - Btt[i] * v[bus_t[i],h] - Btf[i] * c[i,h] - Gtf[i] * s[i,h];


#def i2 var, eq (10), only for goodi2s (c.f. Section 4.5)
			   
subject to i2f_def {i in goodi2, h in 0..T-1}:
        i2f[i,h] = ( ( g[i] * g[i] + b[i] * b[i] ) / ( ratio[i] * ratio[i] ) ) * ( ( v[bus_f[i],h] / (ratio[i] * ratio[i]) ) + v[bus_t[i],h] - ( 2 / ratio[i] ) * ( c[i,h] * cos(phase_angle[i]) + s[i,h] * sin(phase_angle[i]) ) ) + ( b[i] * bshunt[i] / ( ratio[i] * ratio[i] * ratio[i] ) ) * ( ( v[bus_f[i],h] / ratio[i] ) - ( c[i,h] * cos(phase_angle[i]) + s[i,h] * sin(phase_angle[i]) ) ) + ( g[i] * bshunt[i] / ( ratio[i] * ratio[i] * ratio[i] ) ) * ( s[i,h] * cos(phase_angle[i]) - c[i,h] * sin(phase_angle[i]) ) + (bshunt[i] * bshunt[i] * v[bus_f[i],h] / ( 4 * ratio[i] * ratio[i] * ratio[i] * ratio[i] ) );

#badi2 constraints, eqs (19a), (19b)
subject to uppi2 {i in badi2, h in 0..T-1}:
    v[bus_f[i],h] + beta_badi2[i] * v[bus_t[i],h] + gamma_badi2[i] * c[i,h] + zeta_badi2[i] * s[i,h] <= rhs_badi2[i];

subject to lowi2 {i in badi2, h in 0..T-1}:
    v[bus_f[i],h] + beta_badi2[i] * v[bus_t[i],h] + gamma_badi2[i] * c[i,h] + zeta_badi2[i] * s[i,h] >= 0;


#power balance, eqs (1b), (1c)

subject to Pbalance {i in buses, h in 0..T-1}:
       (if i in bus_Gs then (Gs[i] * v[i,h] ) else 0  ) + (sum {j in branches_f[i]} Pf[j,h] ) + (sum {j in branches_t[i]} Pt[j,h] ) = (if card(bus_gens[i]) > 0 then (sum {k in bus_gens[i]} Pg[k,h] ) else 0) - Pd[i,h];

subject to Qbalance {i in buses, h in 0..T-1}:
       (if i in bus_Bs then (- Bs[i] * v[i,h] ) else 0  ) + (sum {j in branches_f[i]} Qf[j,h] ) + (sum {j in branches_t[i]} Qt[j,h] ) = (if card(bus_gens[i]) > 0 then (sum {k in bus_gens[i]} Qg[k,h] ) else 0) - Qd[i];


#jabrs, eq (5)
subject to jabr {i in badi2, h in 0..T-1}: c[i,h] ^2 + s[i,h] ^2 <= v[bus_f[i],h] * v[bus_t[i],h];

#i2s, eq (9)
subject to i2 {i in goodi2, h in 0..T-1}: Pf[i,h] ^2 + Qf[i,h] ^2 <= v[bus_f[i],h] * i2f[i,h];


#limits, eq (1k)

subject to limits_f {i in branches, h in 0..T-1}: Pf[i,h] ^2 + Qf[i,h] ^2 <= U[i] ^2;

subject to limits_t {i in branches, h in 0..T-1}: Pt[i,h] ^2 + Qt[i,h] ^2 <= U[i] ^2;


#ramping constraints, eqs (2a), (2b)

subject to ramp_up {i in gens, h in 0..T-2}: Pg[i,h+1] - Pg[i,h] - rampru[i,h+1] * absPg[i,h] <= 0;

subject to ramp_down {i in gens, h in 0..T-2}: Pg[i,h] - ramprd[i,h+1] * absPg[i,h] - Pg[i,h+1] <= 0;

#abs gens, def absolute value of Pg

subject to abs_gens_plus {i in gens, h in 0..T-1}: Pg[i,h] - absPg[i,h] <= 0;

subject to abs_gens_minus {i in gens, h in 0..T-1}: - Pg[i,h] - absPg[i,h] <= 0;


#def scaledPg, change of vars for Pg of gens with quadcost

subject to scaledPg_def {i in quadgens, h in 0..T-1}: scaledPg[i,h] = alfa[i] * Pg[i,h] - beta[i];





