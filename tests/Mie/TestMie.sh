#!/bin/bash

###################################################################
# check that scuff-scatter is present
###################################################################
CODE=scuff-scatter
which ${CODE} > /dev/null 2>&1
if [ $? -ne 0 ]
then
  echo "could not find ${CODE} executable (check PATH?)"
  exit 1
fi

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
# try to set up timing measurement
##################################################
TIMECMD=
if [ -x `which time` ]
then
  TIMECMD="`which time` -v -o ${RUNTIME}"
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
# if we were invoked with --reference, just create
# and store output files in the `reference` subdirectory,
# then exit
##################################################
FIELDS=EPFile.total
EMTPFT=EMTPFT
DSIPFT=DSIPFT
MOMENTS=Moments
RUNTIME=Runtime
LOG=scuff-scatter.log

if [ ${REFERENCE} -eq 1 ]
then
  mkdir -p reference
  for FILE in ${FIELDS} ${EMTPFT} ${DSIPFT} ${MOMENTS} ${RUNTIME}
  do 
    /bin/mv ${GEOM}.${FILE} reference
    echo "Wrote reference file: reference/${GEOM}.${FILE}"
  done
  /bin/mv ${LOG} reference/${GEOM}.log
  exit 0
fi

##################################################
# check results
##################################################

# if the location of CheckSCUFFData wasn't specified,
# assume it lives one level up in the directory hierarchy
if [ "x${CHECKSCUFFDATA}" == "x" ]
then
  export CHECKSCUFFDATA=../CheckSCUFFData
fi

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

  WALL=`cat ${RUNTIME} | grep 'Elapsed' | cut -f5- -d':'`
  USER=`cat ${RUNTIME} | grep 'User time' | cut -f2- -d':'`
  MEM=`cat ${RUNTIME}  | grep 'Maximum resident set size' | cut -f2- -d':'`

  WALLREF=`cat reference/${TIMEREF} | grep 'Elapsed' | cut -f5- -d':'`
  USERREF=`cat reference/${TIMEREF} | grep 'User time' | cut -f2- -d':'`
  MEMREF=`cat reference/${TIMEREF}  | grep 'Maximum resident set size' | cut -f2- -d':'`

  echo "Wall time: ${WALL} (ref: ${WALLREF})"
  echo "User time: ${USER} (ref: ${USERREF})"
  echo "Memory:    ${MEM}  (ref: ${MEMREF})"
fi

##################################################
##################################################
##################################################
exit ${STATUS}
