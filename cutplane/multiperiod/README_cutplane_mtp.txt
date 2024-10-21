--------

Quick guide on how to run the multi-period cutting-plane script.


1) Edit configuration file 'cutplane_mtp.conf' in directory 'runs':
   Provide the casefilename of the instance;
   Choose the number of time periods T: 4, 12, or 24;
   Remove the string 'addcuts' if running the non warm-started algorithm;
   

   Please check the script 'main_mtp.py' in '/src/' for
   other configuration options (e.g. max number of rounds, time limit,
   tolerances, cut management, type of cuts, write solutions to a .txt
   file, write the model to a .lp file, etc).
   Default values for model and solver parameters can be checked there.
   

2) Navigate to the 'runs' directory to run the command:
   'python3 ../src/main_mtp.py cutplane_mpt.conf'


Example:

We run a warm-started 4-period version of pglib_opf_case10000_goc__api.m
with the following parameters: a time limit of 1200 seconds;
max 20 rounds; 1e-8 feasibility and optimality tolerances;
max 5 iterations without sufficient objective improvement; a
threshold for objective relative improvement of 1e-5; and outer-envelope
cuts Jabr, i2, and limit.

1) The contents of the config file 'cutplane_mtp.conf' are as follows:

   casefilename ../data/pglib_opf_case10000_goc__api.m
   T 4
   addcuts
   max_time 1200
   max_rounds 20
   feastol 1e-8
   optol 1e-8
   ftol_iterates 5
   ftol 1e-5
   jabrcuts
   dropjabr
   i2cuts
   dropi2
   limitcuts
   droplimit   
   END
   
2) 'python3 ../src/main_mtp.py cutplane_mpt.conf'


--------
