#!/bin/bash

. setupSCUFF.sh

#################################################################
#################################################################
#################################################################
L=8
for N in Coarse Fine
do
  ARGS=""
  ARGS="${ARGS} --geometry ${GEODIR}/EFAntenna_L${L}_${N}.scuffgeo"
  ARGS="${ARGS} --PortFile ${PORTDIR}/EFAntenna.ports"
  ARGS="${ARGS} --minFreq  0.5"
  ARGS="${ARGS} --maxFreq  20"
  ARGS="${ARGS} --numFreqs 40"
  ARGS="${ARGS} --ZParameters"
  ARGS="${ARGS} --Eps 2.2 --h 0.794"
  
  ${SCUFF_BIN}scuff-microstrip ${ARGS}
done  
