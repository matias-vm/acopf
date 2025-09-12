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

class OPFException(Exception):
    "DCOPF exception class"

    def _get_message(self):
        return self._message

    def _set_message(self, message):
        self._message = message

    message = property(_get_message, _set_message)

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message
