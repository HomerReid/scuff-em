#!/bin/bash

#################################################################
# SCUFF_BIN and SCUFF_SHARE should be set to the ~bin and ~share directories 
# of the SCUFF-EM installation tree
# FILEBASE should be set to the ~share/scuff-em/examples/MicrostripDevices
#  directory of the SCUFF-EM installation tree; the following attempts
#  to determine its location automatically, but feel free to override
#  by setting BASE appropriately for your system.
#################################################################
SCUFF_BIN=`pkg-config scuff-em --variable=exec-prefix`
if [ "x${SCUFF_BIN}" != "x" ]
then
  SCUFF_BIN=${SCUFF_BIN}/bin
fi

SCUFF_SHARE=`pkg-config scuff-em --variable=datadir`
BASE=${DATADIR}/scuff-em/examples/MicrostripDevices

#################################################################
# set some derived path names
#################################################################
export BASE=${SCUFF_SHARE}/scuff-em/examples/MicrostripDevices
export SCUFF_MESH_PATH=${BASE}/mshFiles
export PORTDIR=${BASE}/portFiles
export GEODIR=${BASE}/scuffgeoFiles

#################################################################
# try to determine the number of processors available and 
# configure OpenMP multithreading accordingly
#################################################################
NPROC=`getconf _NPROCESSORS_ONLN`
if [ "x${NPROC}" != "x" ]
then
  export OMP_NUM_THREADS=${NPROC}
  export GOMP_CPU_AFFINITY="0-$((NPROC-1))"
fi
