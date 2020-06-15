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
# check that scuff-scatter is present
###################################################################
CODE=scuff-scatter
export SCUFF_LOGLEVEL=VERBOSE2
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
GEOM=LossySphere_327
ARGS=""
ARGS="${ARGS} --geometry       ${GEOM}.scuffgeo"
ARGS="${ARGS} --OmegaFile      OmegaFile"
ARGS="${ARGS} --PWDirection    0 0 1"
ARGS="${ARGS} --PWPolarization 0.70710678 0.70710678i 0"
ARGS="${ARGS} --EPFile         EPFile"
ARGS="${ARGS} --MomentFile     ${GEOM}.Moments"
ARGS="${ARGS} --EMTPFTFile     ${GEOM}.EMTPFT"
ARGS="${ARGS} --DSIPFTFile     ${GEOM}.DSIPFT"
ARGS="${ARGS} --DSIRadius      1.5"
ARGS="${ARGS} --DSIPoints      302"

##################################################
##################################################
##################################################
RUNTIME=${GEOM}.Runtime
TIMECMD=`which gtime`
if [ $? -eq 0 ]
then
  ${TIMECMD} -o ${RUNTIME} ${CODE} ${ARGS}
else
  echo "No time command found." > ${RUNTIME}
  ${CODE} ${ARGS}
fi

##################################################
# if we were invoked with --reference, just create
# and store output files in the `reference` subdirectory,
# then exit
##################################################
FIELDS=EPFile.total
EMTPFT=EMTPFT
DSIPFT=DSIPFT
MOMENTS=Moments
LOG=scuff-scatter.log

if [ ${REFERENCE} -eq 1 ]
then
  mkdir -p reference
  for FILE in ${FIELDS} ${EMTPFT} ${DSIPFT} ${MOMENTS}
  do 
    /bin/mv ${GEOM}.${FILE} reference
    echo "Wrote reference file: reference/${GEOM}.${FILE}"
  done
  /bin/mv ${RUNTIME} reference/${RUNTIME}
  /bin/mv ${LOG}     reference/${GEOM}.log
  exit 0
fi

##################################################
# check results
##################################################
STATUS=0
for FILE in ${FIELDS} ${EMTPFT} ${DSIPFT} ${MOMENTS}
do
   ${CHECKSCUFFDATA} --data ${GEOM}.${FILE} --reference reference/${GEOM}.${FILE} --checklist Checklist.${FILE}
   LASTSTATUS=$?
   if [ ${LASTSTATUS} -eq 0 ]
   then
     echo "${FILE} test: success"
   else
     echo "${FILE} test: failure(${LASTSTATUS})"
     STATUS=${LASTSTATUS}
   fi
done

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
