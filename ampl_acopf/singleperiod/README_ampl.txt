--------

Quick guide on how to run the AMPL scripts for single-period ACOPF and SOCP
relaxations:

** Single-period ACOPF **

1) Edit configuration file 'ac.conf' in directory 'runs':
   Provide the casefilename of the instance;

   Please check the script 'main.py' in '/src/' for
   other configuration options. Note that the only nonlinear
   solver set to solve ACOPF is Knitro.


2) Navigate to the 'runs' directory to run the command:
   'python3 ../src/main.py ac.conf'


Example:

We run the single-period version of case_ACTIVSg25k.m with
the solver Knitro with a time limit of 2400 seconds:
	
   casefilename ../data/case_ACTIVSg25k.m
   solver knitroampl
   modfile ../modfiles/ac.mod
   max_time 2400
   END


** Single-period SOCP **

1) Edit config file 'socp.conf' in directory 'runs':
   Provide the casefilename of the instance;
   Provide the string "solver <solvername>"

   Please check the script 'main.py' in '/src/' for
   other configuration options.

2) Navigate to the 'runs' directory to run the command:
   'python3 ../src/main.py socp.conf'


Example:

We run the single-period version of case9241pegase.m Jabr SOCP with the
solver Knitro with a time limit of 1000 seconds.

1) The contents of the config file 'socp_mtp.conf' are as follows:

   casefilename ../data/case9241pegase.m
   solver knitroampl
   modfile ../modfiles/jabr.mod
   max_time 1000
   END

2) 'python3 ../src/main.py socp.conf'

**
To change Mosek parameters the file '/opt/mosek/10/mosekopt.par' needs to be edited.
**

**
Please check the modfiles 'ac_mtp_definedvars' and 'i2plus_mosek_mtp.mod for model
descriptions with references to equations in our paper.
**

**
Knitro and Gurobi share the same modfile for the Jabr and i2 SOCPs:
'../modfiles/jabr.mod' and '../modfiles/i2.mod'. 
If Mosek is called, then please use '../modfiles/jabr_mosek.mod' and
'../modfiles/i2_mosek.mod'.
**

**
The data provided corresponds to a subset of the largest instances
from PGLIB, the Pegase project, and ACTIVSg cases. The two type of
single-period perturbed cases are:

    Loads perturbed – Gaussian deltas: suffix '_n_1_1.m'
    Transmission line with largest flow switched off: suffix '_line.m'
**

--------
