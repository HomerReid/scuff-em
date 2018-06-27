###############################################################
# Some general-purpose setup for python-driven SCUFF-EM.
#
# Usage: Just include the line
#
#   execfile('SetupSCUFF.py')
#
# near the top of your python script (or at the command line).
# 
# Optionally, if you set the variable "scuffExample" to the
# name of one of the folders in share/scuff-em/examples
# (i.e. scuffExample="MicrostripDevices") before the 'execfile'
# line then various environment variables will be set to
# make it easier for SCUFF-EM to find input files.
#
# Homer Reid 20180422
###############################################################
from __future__ import print_function;

import os;
import subprocess;
import commands;
import numpy as np;
import time;
import math;
import argparse;

import scuff;

##############################################################
# Try to determine automatically the installed location of
# SCUFF-EM. You can bypass this by presetting the environnment
# variable SCUFF to the "--prefix" you specified when
# installing SCUFF-EM (default is /usr/local.)
##############################################################
# first attempt: look at SCUFF environment variable
if os.getenv("SCUFF") is not None:
    SCUFF=os.environ["SCUFF"]

# next attempt: use pkg-config
if 'SCUFF' not in locals() and 'SCUFF' not in globals():

    os.environ["SCUFF_MESH_PATH"] = DirBase + "mshFiles"

    SCUFF="/usr/local" # default
    try:
        SCUFF=commands.getoutput("pkg-config scuff-em --variable=exec_prefix");
    except:
        pass

if 'SCUFF' in locals() or 'SCUFF' in globals():
    os.environ["SCUFF"]=SCUFF

SCUFF_EXAMPLES=SCUFF + "/share/scuff-em/examples"
##############################################################
# If the user set the variable scuffExample to the name of one
# of the scuff-em example (i.e. scuffExample="MicrostripDevices")
# then set some environment variables so SCUFF-EM knows where to
# look for input files.
##############################################################
if 'scuffExample' in locals() or 'scuffExample' in globals():
    DirBase = SCUFF_EXAMPLES + "/" + scuffExample + "/"
    os.environ["SCUFF_MESH_PATH"] = DirBase + "mshFiles"
    os.environ["SCUFF_GEO_PATH"]  = DirBase + "scuffgeoFiles"
    os.environ["SCUFF_DATA_PATH"] = DirBase + "portFiles"

print("SCUFF_MESH_PATH = ", os.environ["SCUFF_MESH_PATH"])

##############################################################
# try to determine number of CPU cores and ask for that many
# threads when running SCUFF
##############################################################
import multiprocessing;
try:
    NumCores=multiprocessing.cpu_count()
    os.environ["OMP_NUM_THREADS"]="{}".format(NumCores)
    os.environ["GOMP_CPU_AFFINITY"]="{}-{}".format(0,NumCores-1)
    print("Using {} CPU cores.".format(NumCores));
except:
    print("Failed to configure multithreading.");

##############################################################
# do some preliminary SCUFF-EM bookkeeping normally handled   
# at the top of SCUFF-EM API codes                            
##############################################################
scuff.SetLogFileName("pyscuff.log");
scuff.InstallHRSignalHandler();
#os.environ["SCUFF_LOGLEVEL"]="VERBOSE2"
