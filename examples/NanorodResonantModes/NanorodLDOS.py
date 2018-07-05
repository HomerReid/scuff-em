############################################################## 
# python code for LDOS near gold nanorod(s)
##############################################################
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
# user-tweakable parameters
##############################################################
OmegaRange  = np.linspace(5.0,8.0,50)
Res="648"     # resolution [648 | 1860 | 3612]
Dimer=True;  # analyze dimer instead of monomer geometry

##############################################################
# try to figure out where SCUFF-EM example files are installed
############################################################## 
SCUFF_PREFIX="/usr/local"
if os.getenv("SCUFF_PREFIX") is not None:
    SCUFF_PREFIX=os.environ["SCUFF_PREFIX"]
else:
    (status,result)=commands.getstatusoutput("pkg-config scuff-em --variable=prefix");
    if (status==0):
        SCUFF_PREFIX=result

BaseDir=SCUFF_PREFIX + "/share/scuff-em/examples/NanorodResonantModes"
os.environ["SCUFF_MESH_PATH"] = BaseDir + "/mshFiles"
os.environ["SCUFF_MATPROP_PATH"] = BaseDir

###################################################
# create the solver, define the geometry
###################################################
scuff.SetLogFileName("NanorodLDOS.log")
scuff.InstallHRSignalHandler();
Name = ("Nanorods" if Dimer else "Nanorod") + "_" + Res
Solver=scuff.scuffSolver(Name)

MeshFile="Nanorod_" + Res + ".msh";

Solver.SetMediumPermittivity(2.25)

Material="Gold"
Solver.AddObject(MeshFile, Material, "Nanorod")
if Dimer:
    Transformation="DISPLACED 0 0.05 0"
    Solver.AddObject(MeshFile, Material, "Nanorod2", Transformation)

SourcePoint = [0.0, 0.0, 0.09]

###################################################
###################################################
###################################################
LDOSFile=open(Name + ".LDOS", "w")
for Omega in OmegaRange:

    Solver.AssembleSystemMatrix(Omega)

    GEDiag=[0,0,0]
    for Mu in range(0,3):

       # create Mu-directed dipole source
       SourceOrientation=[0,0,0]
       SourceOrientation[Mu]=1
       PS=scuff.PointSource(SourcePoint, SourceOrientation)

       # solve scattering problem
       Solver.Solve(PS)

       # get scattered E, H fields at source location and 
       # accumulate E[Mu] component
       EH=Solver.GetFields(SourcePoint,"scattered")
       GEDiag[Mu]=EH[Mu]

    LDOS=sum(GEDiag)
    LDOSFile.write("%e %e %e " % (Omega,LDOS.real, LDOS.imag))
    LDOSFile.write("%e %e "    % (GEDiag[0].real,GEDiag[0].imag))
    LDOSFile.write("%e %e "    % (GEDiag[1].real,GEDiag[1].imag))
    LDOSFile.write("%e %e\n"   % (GEDiag[2].real,GEDiag[2].imag))
    LDOSFile.flush()

LDOSFile.close()
