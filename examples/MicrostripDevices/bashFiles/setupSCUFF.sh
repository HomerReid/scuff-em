#!/bin/bash

#################################################################
# The following is an attempt to discern automatically the root
# of the SCUFF-EM installation tree on your system; you can
# bypass this by presetting the environment variable SCUFF
# to the `--prefix` you specified when configuring your SCUFF-EM
# build (default is /usr/local/share
#SCUFF_BIN and SCUFF_SHARE should be set to the ~bin and ~share directories of
# of the SCUFF-EM installation tree.
# FILEBASE should be set to the ~share/scuff-em/examples/MicrostripDevices
#  directory of the SCUFF-EM installation tree; the following attempts
#  to determine its location automatically, but feel free to override
#  by setting BASE appropriately for your system.
#################################################################
if [ "x${SCUFF}" == "x" ]
then
  SCUFF=`pkg-config scuff-em --variable=exec-prefix`
  if [ "x${SCUFF}" == "x" ]
  then
    SCUFF=/usr/local # unlikely this will work, since then the pkg-config would have worked
  fi
fi

#################################################################
# set some derived path names
#################################################################
BASE=${SCUFF}/share/scuff-em/examples/MicrostripDevices
export SCUFF_MESH_PATH=${BASE}/mshFiles
export PORTDIR=${BASE}/portFiles
export GEODIR=${BASE}/scuffgeoFiles

SCUFF_BIN=
if [ "x${SCUFF}" != "x" ]
then
  export SCUFF_BIN="${SCUFF}/bin/"
fi

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
