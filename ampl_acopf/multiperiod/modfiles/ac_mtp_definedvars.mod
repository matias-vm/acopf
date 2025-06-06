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

param T; # number of time periods
param fixedcost {i in gens}; #constant term objective
param lincost {i in gens}; #coeff linear terms objective
param quadcost {i in gens}; #coeff quad terms objective

#branch parameters
param Gtt {i in branches}; 
param Btt {i in branches};
param Gff {i in branches};
param Bff {i in branches};
param Gtf {i in branches};
param Btf {i in branches};
param Gft {i in branches};
param Bft {i in branches};

param bus_f {i in branches}; #from buses of each branch
param bus_t {i in branches}; #to buses of each branch
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
param theta_max {i in buses}; #max angle
param theta_min {i in buses}; #min angle
param maxangle {i in branches}; #max angle diff
param minangle {i in branches}; #min angle diff

param Vinit {i in buses, h in 0..T-1}; #flat-start, initial point
param thetadiffinit {i in branches, h in 0..T-1}; #flat-start, initial point

param rampru {i in gens, h in 0..T-1}; #ramp up rates
param ramprd {i in gens, h in 0..T-1}; #ramp down rates

#VARIABLES
var v {i in buses, h in 0..T-1} >= Vmin[i], <= Vmax[i], := Vinit[i,h]; # eq (1j)
var theta {i in buses, h in 0..T-1} >= theta_min[i], <= theta_max[i];
var thetadiff {i in branches, h in 0..T-1} >= minangle[i], <= maxangle[i], := thetadiffinit[i,h];
var Pg {i in gens, h in 0..T-1} >= Pmin[i], <= Pmax[i];
var Qg {i in gens, h in 0..T-1} >= Qmin[i], <= Qmax[i];
var absPg {i in gens, h in 0..T-1} >= 0;


#DEFINED VARS

#defP, active power flows eqs (1d), (1e)

var Pf{i in branches, h in 0..T-1} =
                   Gff[i] * v[bus_f[i],h] * v[bus_f[i],h] + Gft[i] * v[bus_f[i],h] * v[bus_t[i],h] * cos(thetadiff[i,h]) + Bft[i] * v[bus_f[i],h] * v[bus_t[i],h] * sin(thetadiff[i,h]);

var Pt{i in branches, h in 0..T-1} =
                   Gtt[i] * v[bus_t[i],h] * v[bus_t[i],h] + Gtf[i] * v[bus_f[i],h] * v[bus_t[i],h] * cos(thetadiff[i,h]) - Btf[i] * v[bus_f[i],h] * v[bus_t[i],h] * sin(thetadiff[i,h]);

#defQ, reactive power flows eqs (1f), (1g)

var Qf{i in branches, h in 0..T-1} =
                   - Bff[i] * v[bus_f[i],h] * v[bus_f[i],h] - Bft[i] * v[bus_f[i],h] * v[bus_t[i],h] * cos(thetadiff[i,h]) + Gft[i] * v[bus_f[i],h] * v[bus_t[i],h] * sin(thetadiff[i,h]);

var Qt{i in branches, h in 0..T-1} =
                   - Btt[i] * v[bus_t[i],h] * v[bus_t[i],h] - Btf[i] * v[bus_f[i],h] * v[bus_t[i],h] * cos(thetadiff[i,h]) - Gtf[i] * v[bus_f[i],h] * v[bus_t[i],h] * sin(thetadiff[i,h]);

#OBJECTIVE, eq (1a)

minimize total_cost: sum {i in gens, h in 0..T-1} (fixedcost[i] + lincost[i] * Pg[i,h] + quadcost[i] * Pg[i,h] ^2);

#power balance, eqs (1b) and (1c)

subject to Pbalance {i in buses, h in 0..T-1}:
       (if i in bus_Gs then (Gs[i] * v[i,h] * v[i,h] ) else 0  ) + (sum {j in branches_f[i]} Pf[j,h] ) + (sum {j in branches_t[i]} Pt[j,h] ) = (if card(bus_gens[i]) > 0 then (sum {k in bus_gens[i]} Pg[k,h] ) else 0) - Pd[i,h];

subject to Qbalance {i in buses, h in 0..T-1}:
        (if i in bus_Bs then (- Bs[i] * v[i,h] * v[i,h] ) else 0  ) + (sum {j in branches_f[i]} Qf[j,h] ) + (sum {j in branches_t[i]} Qt[j,h] ) = (if card(bus_gens[i]) > 0 then (sum {k in bus_gens[i]} Qg[k,h] ) else 0) - Qd[i];


#limits, capacity constraints, eq (1k)

subject to limits_f {i in branches, h in 0..T-1}: Pf[i,h] ^2 + Qf[i,h] ^2 <= U[i] ^2;

subject to limits_t {i in branches, h in 0..T-1}: Pt[i,h] ^2 + Qt[i,h] ^2 <= U[i] ^2;

#angle diff

subject to thetadiff_def {i in branches, h in 0..T-1}: thetadiff[i,h] = theta[bus_f[i],h] - theta[bus_t[i],h];

#ramping constraints, eqs (2a), (2b)

subject to ramp_up {i in gens, h in 0..T-2}: Pg[i,h+1] - Pg[i,h] - rampru[i,h+1] * absPg[i,h] <= 0;

subject to ramp_down {i in gens, h in 0..T-2}: Pg[i,h] - ramprd[i,h+1] * absPg[i,h] - Pg[i,h+1] <= 0;

#abs gens

subject to abs_gens_plus {i in gens, h in 0..T-1}: Pg[i,h] - absPg[i,h] <= 0;

subject to abs_gens_minus {i in gens, h in 0..T-1}: - Pg[i,h] - absPg[i,h] <= 0;
