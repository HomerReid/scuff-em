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
# set up arguments for SCUFF-EM run and subsequent result checking
###################################################################
CODE=scuff-transmission
ARGS=""
ARGS="${ARGS} --geometry    E10HalfSpace_40.scuffgeo"
ARGS="${ARGS} --Omega       0.1"
ARGS="${ARGS} --Omega       1.0"
ARGS="${ARGS} --ThetaMin    0"
ARGS="${ARGS} --ThetaMax    85"
ARGS="${ARGS} --ThetaPoints 19"
ARGS="${ARGS} --ZAbove      1.0"

CHECKLIST=Fresnel.Checklist
DATAFILE=E10HalfSpace_40.transmission
DATAREF=${DATAFILE}.reference
TIMEFILE=E10HalfSpace_40.timing
TIMEREF=${TIMEFILE}.reference

##################################################
# try to set up timing measurement
##################################################
TIMECMD=
if [ -x `which time` ]
then
  TIMECMD="`which time` -v -o ${TIMEFILE}"
fi

# if the location of CheckSCUFFData wasn't specified,
# assume it lives one level up in the directory hierarchy
if [ "x${CHECKSCUFFDATA}" == "x" ]
then
  export CHECKSCUFFDATA=../CheckSCUFFData
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

##################################################
# run the code
##################################################
${TIMECMD} ${CODE} ${ARGS}

##################################################
# if we were invoked with --reference, just store
# output files as reference files and exit       
##################################################
if [ ${REFERENCE} -eq 1 ]
then
  /bin/mv ${DATAFILE} ${DATAREF}
  /bin/mv ${TIMEFILE} ${TIMEREF}
  echo "Wrote reference files ${DATAREF} and ${TIMEREF}."
  exit 0
fi

##################################################
# check results
##################################################
${CHECKSCUFFDATA} --data ${DATAFILE} --reference ${DATAREF} --checklist ${CHECKLIST}
Status=$?

##################################################
# extract timing info if available 
##################################################
if [ -r ${TIMEFILE} ]
then

  WALL=`cat ${TIMEFILE} | grep 'Elapsed' | cut -f5- -d':'`
  USER=`cat ${TIMEFILE} | grep 'User time' | cut -f2- -d':'`
  MEM=`cat ${TIMEFILE}  | grep 'Maximum resident set size' | cut -f2- -d':'`

  WALLREF=`cat ${TIMEREF} | grep 'Elapsed' | cut -f5- -d':'`
  USERREF=`cat ${TIMEREF} | grep 'User time' | cut -f2- -d':'`
  MEMREF=`cat ${TIMEREF}  | grep 'Maximum resident set size' | cut -f2- -d':'`

  echo "Wall time: ${WALL} (ref: ${WALLREF})"
  echo "User time: ${USER} (ref: ${USERREF})"
  echo "Memory:    ${MEM}  (ref: ${MEMREF})"
fi

##################################################
##################################################
##################################################
exit ${STATUS}
