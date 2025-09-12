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

def break_exit(foo):

    stuff = input("("+foo+") break> ")
    if stuff == 'x' or stuff == 'q':
        sys.exit("bye")
