################################################################################                            
##                                                                            ##                            
## This code was adapted and extended from the Gurobi OptiMod                 ##                            
## 'Optimal Power Flow', written by Daniel Bienstock. It is being maintained  ##                            
## by Matias Villagra.                                                        ##                            
##                                                                            ##                            
## Please report any bugs or issues to                                        ##                            
##                        mjv2153@columbia.edu                                ##                            
##                                                                            ##                            
## Jul 2024                                                                   ##                           
################################################################################

import sys
import time
from socket import gethostname
from opfexception import OPFException

class Logger:
    """Class to handle logging of the OPF solve beyond Gurobi output"""

    def __init__(self, logfilename):
        self.logfile = open(logfilename,"w")
        self.screen  = 1
        self.log     = 1
        localtime    = time.asctime(time.localtime(time.time()))
        self.joint("Starting log: %s\n" % localtime)
        self.joint("Running on: " + gethostname() + "\n")
        self.printversion()

    # Print date and close logfile
    def close_log(self):
        localtime = time.asctime(time.localtime(time.time()))
        self.joint("Closing log: %s\n" % localtime)
        self.logfile.close()

    # Print output to screen and/or logfile
    def joint(self, message, *args):
        write = 1
        for arg in args:
            write = arg

        # Write to logfile
        if self.log and write == 1:
            self.logfile.write(message)
            self.logfile.flush()

        # Print to screen
        if self.screen and write == 1:
            print(message, end = ''),
            sys.stdout.flush() 

    def both_on(self):
        self.screen = self.log = 1

    def both_off(self):
        self.screen = self.log = 0

    def screen_on(self):
        self.screen = 1

    def screen_off(self):
        self.screen = 0

    def log_on(self):
        self.log = 1

    def log_off(self):
        self.log = 0

    # Show error message and raise exception
    def raise_exception(self, message):
        self.log    = 1
        self.screen = 0
        if message:
            self.joint("\n" + message)
        self.close_log()
        raise OPFException("\n DCOPFException: " + message + "Encountered an error. Quitting")


    def printversion(self):
        self.joint("Version Sun.Sep..08.201000.2024@mpc-linux-01\n\n")


