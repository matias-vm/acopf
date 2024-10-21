--------

Quick guide on how to run the single-period cutting-plane script.


1) Edit configuration file 'cutplane.conf' in directory 'runs':
   Provide the casefilename of the instance;
   Remove the string 'addcuts' if running the non warm-started algorithm;
   
   Please check the script 'main.py' in '/src/' for
   other configuration options (e.g. max number of rounds, time limit,
   tolerances, cut management, type of cuts, write solutions to a .txt
   file, write the model to a .lp file, etc).
   Default values for model and solver parameters can be checked there.


2) Navigate to the 'runs' directory to run the command:
   'python3 ../src/main.py cutplane.conf'


Example:

We run the warm-started single-period instance of case_ACTIVSg10k.m
with the following parameters: a time limit of 1000 seconds;
max 100 rounds; max 5 iterations without sufficient objective improvement;
and a threshold for objective relative improvement of 1e-5; and outer-envelope
Jabr, i2 and limit cuts.

1) The contents of the config file 'cutplane_mtp.conf' are as follows:

   casefilename ../data/case_ACTIVSg10k.m
   addcuts
   max_time 1000
   max_rounds 20
   ftol_iterates 5
   ftol 1e-5
   jabrcuts
   dropjabr
   i2cuts
   dropi2
   limitcuts
   droplimit
   END

2) 'python3 ../src/main.py cutplane.conf'


***
For a detailed description of the model variables and constraints please                                
check 'cutplane_mtp.py' in '../../multiperiod/src/', which corresponds
to the multi-period version of this model.

The data provided corresponds to a subset of the largest instances
from PGLIB, the Pegase project, and ACTIVSg cases. The two type of
single-period perturbed cases are:

    Loads perturbed – Gaussian deltas: suffix '_n_1_1.m'
    Transmission line with largest flow switched off: suffix '_line.m'
    
***

--------
