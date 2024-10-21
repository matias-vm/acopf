###############################################################################
##                                                                           ##
## This code was written by Daniel Bienstock and is being maintained by      ##
## Matias Villagra, PhD Student in Operations Research @ Columbia,           ##
## supervised by Daniel Bienstock.                                           ##
##                                                                           ##           
## For code readability, we make references to equations in:                 ##
## [1] D. Bienstock, and M. Villagra, Accurate Linear Cutting-Plane          ##
## Relaxations for ACOPF, arXiv:2312.04251v2, 2024                           ##
##                                                                           ##
## Please report any bugs or issues to: mjv2153@columbia.edu                 ## 
##                                                                           ##
## Jul 2024                                                                  ##
###############################################################################

# -*- coding: utf-8 -*-
"""

""" 
import sys
import time

from log import danoLogger

def stateversion(log):
    log.joint("Version <Wed.Jul.10.130144.2024@cool1>\n")

