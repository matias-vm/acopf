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
from myutils import breakexit
from socket import gethostname

#this class handles logging
class danoLogger:
  def __init__(self, logfilename):
    self.logfile = open(logfilename,"a+")
    self.screen = 1
    self.log = 1
    localtime= time.asctime(time.localtime(time.time()))
    self.joint("starting log: %s\n" % localtime)
    self.joint("running on: " + gethostname() + "\n\n")
   


  def closelog(self):
    localtime= time.asctime(time.localtime(time.time()))
    self.logfile.write("\nclosing log: %s\n" % localtime)
    self.logfile.close()

  def joint(self, mystring, *args):
    write = 1
    for arg in args:
        write = arg
    if self.log and write==1:
      self.logfile.write(mystring)
      self.logfile.flush()
    if self.screen and write==1:
      print(mystring, end=''),
      sys. stdout. flush() 
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

  def stateandquit(self, stuff):
    self.log = self.screen = 1
    self.joint(stuff + "\n")
    self.closelog()
    sys.exit("\nquitting\n")    
     




