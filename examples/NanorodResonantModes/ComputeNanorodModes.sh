#!/bin/bash

##################################################
# Simple bash script that uses SCUFF-SPECTRUM to
# compute resonant frequencies and plot field
# patterns for resonant modes of gold nanorod monomer/dimer
#
# Command-line options:
#
#  --Omega0 6.5-0.3i             (center of contour)
#  --RX 0.1                      (horizontal radius of contour ellipse)
#  --RY 0.1                      (vertical radius of contour ellipse)
#  --NQ 10                       (number of quadrature points)
#  --NL 5                        (estimated number of eigenvalues)
#  --Resolution [648|1860|3612]  (fineness of nanorod mesh)
#  --Dimer                       (use nanorod dimer instead of monomer)
##################################################
OMEGA0=6.5-0.3i
RX=0.25
RY=0.1
NQ=20
NL=5
RES=648
DIMER=0
while [ $# -ne 0 ]
do
  case "${1}" in
    "--Omega0")      OMEGA0=$2; shift 2; ;;
    "--RX")          RX=$2;     shift 2; ;;
    "--RY")          RY=$2;     shift 2; ;;
    "--NQ")          NQ=$2;     shift 2; ;;
    "--NL")          NL=$2;     shift 2; ;;
    "--Resolution")  RES=$2;    shift 2; ;;
    "--Dimer")       DIMER=1;   shift 1; ;;
    *) echo "unknown argument $1"; exit 1; ;;
  esac
done

##################################################
# If you installed SCUFF-EM somewhere other than
# the default destination (/usr/local/), set the
# SCUFF_PREFIX environment variable to the installation
# path (i.e. the `--prefix`) you specified.
##################################################
export SCUFF_PREFIX=${SCUFF_PREFIX:-/usr/local}

export SCUFF_BEYN_VERBOSE=1
export SCUFF_BEYN_RANK_TOL=1.0e-4
export SCUFF_BEYN_RES_TOL=1.0e-4

##################################################
# setup some SCUFF-EM directory paths
##################################################
export BASE=${SCUFF_PREFIX}/share/scuff-em/examples/NanorodResonantModes
export SCUFF_MESH_PATH=${BASE}/mshFiles
export SCUFF_GEO_PATH=${BASE}/scuffgeoFiles

##################################################
# try to use all available CPU cores
##################################################
CORES=`nproc`
if [ $? -eq 0 ]
then
  export OMP_NUM_THREADS="${CORES}"
  export GOMP_CPU_AFFINITY="0-$((CORES-1))"
fi

##################################################
# set some SCUFF-EM environment variables
##################################################
export SCUFF_MATRIX_2018=1

##################################################
##################################################
##################################################
if [ "x${DIMER}" == "x1" ]
then
  GEOM=GoldNanorods_${RES}
  LY=0.120 
  NY=24
else
  GEOM=GoldNanorod_${RES}
  LY=0.060 
  NY=12
fi

ARGS=""
ARGS="${ARGS} --geometry ${GEOM}.scuffgeo"
ARGS="${ARGS} --omega0 ${OMEGA0}"
ARGS="${ARGS} --RX ${RX}"
ARGS="${ARGS} --RY ${RY}"
ARGS="${ARGS} --N ${NQ}"
ARGS="${ARGS} --L ${NL}"
ARGS="${ARGS} --FVScreen  0.03 -0.03 -0.10 0    ${LY} 0  0 0 0.2 ${NY} 40"
ARGS="${ARGS} --FVScreen -0.03  0.03 -0.10 0.06 0     0  0 0 0.2 12 40"
ARGS="${ARGS} --PlotContours"
ARGS="${ARGS} --PlotSurfaceCurrents"
ARGS="${ARGS} --PlotSurfaceFields"
${SCUFF_PREFIX}/bin/scuff-spectrum ${ARGS}
