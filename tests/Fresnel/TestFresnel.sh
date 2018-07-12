#!/bin/bash

###################################################################
# process arguments ###############################################
###################################################################
REFERENCE=0
if [ $# == 1 ]
then
  if [ "x${1}" == "x--reference" ]
  then
    REFERENCE=1
    echo "Reference run to generate data files."
  else
    echo "unknown command-line argument ${1} (aborting)"
    exit 1
  fi
fi

###################################################################
# check that scuff-transmission is present
###################################################################
CODE=scuff-transmission
which ${CODE} > /dev/null 2>&1
if [ $? -ne 0 ]
then
  echo "could not find ${CODE} executable (check PATH?)"
  exit 1
fi

###################################################################
# if no explicit location for CheckSCUFFData was furnished, look
# for it in various places
###################################################################
if [ "x${CHECKSCUFFDATA}" == "x" ]
then
  DIR=`which CheckSCUFFData`
  if [ $? -eq 0 ]
  then
    export CHECKSCUFFDATA=${DIR}
  fi
fi

if [ "x${CHECKSCUFFDATA}" == "x" ]
then
  SCUFFDATADIR=`pkg-config scuff-em --variable=datadir`
  if [ $? -eq 0 ]
  then
    export PATH=${PATH}:${SCUFFDATADIR}/tests
  fi

  if [ "x${CHECKSCUFFDATA}" == "x" ]
  then
    echo "could not find CheckSCUFFData executable (set CHECKSCUFFDATA environment variable)"
    exit 1
  fi
fi

##################################################
# try to set the number of CPU cores as high as possible
# (if it is not already set)
##################################################
if [ "x${OMP_NUM_THREADS}" == "x" ]
then 
  CORES=`getconf _NPROCESSORS_ONLN`
  if [ "x${CORES}" != "x" ]
  then
    export OMP_NUM_THREADS=${CORES}
    export GOMP_CPU_AFFINITY=0-$((CORES-1))
    echo "Using ${CORES} CPU cores (${GOMP_CPU_AFFINITY})"
  fi
fi

###################################################################
# set up arguments for SCUFF-EM run and subsequent result checking
###################################################################
export SCUFF_LOGLEVEL="VERBOSE2"
GEOM="E10HalfSpace_40"
ARGS=""
ARGS="${ARGS} --geometry    ${GEOM}.scuffgeo"
ARGS="${ARGS} --Omega       0.1"
ARGS="${ARGS} --Omega       1.0"
ARGS="${ARGS} --ThetaMin    0"
ARGS="${ARGS} --ThetaMax    85"
ARGS="${ARGS} --ThetaPoints 19"
ARGS="${ARGS} --ZAbove      1.0"

##################################################
##################################################
##################################################
TRANSMISSION=${GEOM}.transmission
RUNTIME=${GEOM}.Runtime

TIMECMD=`which time`
if [ $? -eq 0 ]
then
  ${TIMECMD} -o ${RUNTIME} ${CODE} ${ARGS}
else
  echo "No time command found." > ${RUNTIME}
  ${CODE} ${ARGS}
fi

##################################################
# if we were invoked with --reference, just store
# output files as reference files and exit       
##################################################
if [ ${REFERENCE} -eq 1 ]
then
  mkdir -p reference
  /bin/mv ${TRANSMISSION} reference/${TRANSMISSION}
  /bin/mv ${RUNTIME} reference/${RUNTIME}
  echo "Wrote reference files reference/${TRANSMISSION} and reference/${RUNTIME}."
  exit 0
fi

##################################################
# check results
##################################################
CHECKLIST=Fresnel.Checklist
${CHECKSCUFFDATA} --data ${TRANSMISSION} --reference reference/${TRANSMISSION} --checklist ${CHECKLIST}
STATUS=$?

##################################################
# extract timing info if available 
##################################################
if [ -r ${RUNTIME} ]
then
  echo "Timing (this run): "
  cat ${RUNTIME}
  if [ -r reference/${RUNTIME} ]
  then
    echo "Timing (reference): "
    cat reference/${RUNTIME}
  fi
fi

##################################################
##################################################
##################################################
exit ${STATUS}
