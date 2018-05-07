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
import scuff;

##############################################################
# Try to determine automatically the installed location of
# the SCUFF-EM examples and set it as SCUFF_EXAMPLES.
# You can bypass this by presetting SCUFF_EXAMPLES to be
# to ${prefix}/share/scuff-em/examples where ${prefix} is the
# installation prefix you specified when installing SCUFF-EM
# (default is /usr/local.)
##############################################################
# first attempt: use pkg-config
if 'SCUFF_EXAMPLES' not in locals() and 'SCUFF_EXAMPLES' not in globals():

    SCUFF_EXAMPLES="/usr/local/share/scuff-em/examples" # default
    try:
        SCUFF_DATA=commands.getoutput("pkg-config scuff-em --variable=datadir");
        if len(SCUFF_DATA)!=0:
            SCUFF_EXAMPLES=SCUFF_DATA + "/scuff-em/examples"
    except:
        pass

if 'SCUFF_EXAMPLES' in locals() or 'SCUFF_EXAMPLES' in globals():
    os.environ["SCUFF_EXAMPLES"]=SCUFF_EXAMPLES

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
os.environ["SCUFF_LOGLEVEL"]="VERBOSE2"
