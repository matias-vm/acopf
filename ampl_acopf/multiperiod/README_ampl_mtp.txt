--------

Quick guide on how to run the AMPL scripts for multi-period ACOPF and SOCP
relaxations:

** Multi-period ACOPF **

1) Edit configuration file 'ac_mtp.conf' in directory 'runs':
   Provide the casefilename of the instance;
   Choose the number of time periods T: 4, 12, or 24;
   The script runs heuristic3 by default (see below):
	If the string "heuristic1" is given, the script runs heuristic1,
      	If the string "heuristic2" is given, then it runs heuristic2.

   Some heuristics for finding AC primal bounds using Knitro:

   heuristic1: We give to Knitro the full Multi-Time Period Formulation
   heuristic2: We give to Knitro one period at a time and impose
   	       the ramping constraints using generation from the previous
               time period
   heuristic3: First we run heuristic1. If it fails, relax tolerances and
               run heuristic2

   Please check the script 'main_mtp.py' in '/src/' for
   other configuration options. Note that the only nonlinear
   solver set to solve ACOPF is Knitro.


2) Navigate to the 'runs' directory to run the command:
   'python3 ../src/main_mtp.py ac_mtp.conf'


Example:

We run a 4-period version of case_ACTIVSg10k.m with heuristic1
the solver Knitro with a time limit of 2400 seconds:
	
   casefilename ../data/case_ACTIVSg10k.m
   T 4
   heuristic1
   max_time 2400
   END


** Multi-period SOCP **

1) Edit config file 'socp_mtp.conf' in directory 'runs':
   Provide the casefilename of the instance;
   Provide the string 'jabr' to run the Jabr SOCP, while
   the string 'i2plus' runs the i2plus SOCP relaxation;
   Choose the number of time periods T: 4, 12, or 24;
   Provide the string "solver <solvername>" (gurobi, mosek, or knitroampl)

   Please check the script 'main_mtp.py' in '/src/' for
   other configuration options.

2) Navigate to the 'runs' directory to run the command:
   'python3 ../src/main_mtp.py socp_mpt.conf'


Example:

We run a 4-period version of case9241pegase.m i2+ SOCP with the
solver Gurobi with a time limit of 2400 seconds.

1) The contents of the config file 'socp_mtp.conf' are as follows:

   casefilename ../data/case9241pegase.m
   T 4   
   solver gurobi
   i2plus
   max_time 2400
   END

2) 'python3 ../src/main_mtp.py socp_mpt.conf'

**
Please check the modfiles 'ac_mtp_definedvars.mod' and 'i2plus_mosek_mtp.mod' for model
descriptions with references to equations in our paper.
**

**
To change Mosek parameters the file '/opt/mosek/10/mosekopt.par' needs to be edited.
**

**
The data provided corresponds to 4, 12, and 24 time periods (c.f. directories
/data/ramprates and /data/mtploads) and a subset of the largest instances
from PGLIB, the Pegase project, and ACTIVSg cases.
**

--------
