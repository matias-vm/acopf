--------

Quick guide on how to run the multi-period DCOPF script.


1) Edit configuration file 'dcopf_mtp.conf' in directory 'runs':
   Provide the casefilename of the instance;
   Choose the number of time periods: 4, 12, or 24.   
   Feasibility and optimality tolerances are hardcoded (1e-8)

2) Navigate to the 'runs' directory to run the command:
   'python3 ../src/main_mtp.py dcopf_mtp.conf'


Example:

We run a 12-period DCOPF instance case1354pegase.m

1) The contents of the config file 'dcopf_mtp.conf' are as follows:

   casefilename ../data/case1354pegase.m
   T 12
   END

2) 'python3 ../src/main_mtp.py dcopf_mtp.conf'


--------
