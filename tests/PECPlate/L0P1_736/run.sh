#!/bin/bash

#SBATCH -c 20
#SBATCH --partition=hns
#SBATCH --time=47:59:59

export SCUFF_MATRIX_2018=1

scuff-caspol --geometry AuPlate_736.scuffgeo --EPFile ../EPFile --Atom Rb --XiFile ../XiFile --BZIMethod Polar --BZIOrder 3111 --BZSymmetryFactor 8
