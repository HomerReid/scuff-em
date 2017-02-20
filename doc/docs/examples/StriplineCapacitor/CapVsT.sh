#!/bin/bash

CODE=scuff-static
GEOM=StriplineCapacitor

#################################################################
#################################################################
#################################################################
N=4   # meshing fineness
W=0.5 # trace width
/bin/rm ${GEOM}_W${W}_N${N}.CapVsT
for T in `cat TFile`
do
  #################################################################
  # generate geometry mesh for this thickness                      
  #################################################################
  gmsh -2 -setnumber N ${N} -setnumber T ${T} -setnumber W ${W} ${GEOM}.geo
  NUMEDGES=`scuff-analyze --mesh ${GEOM}.msh  | grep 'interior edges' | head -1 | cut -f2 -d' '`
  echo "T=${T}: NumEdges=${NUMEDGES}"

  #################################################################
  # get capacitance matrix for this geometry ######################
  #################################################################
  scuff-static --geometry ${GEOM}.scuffgeo --CapFile ${GEOM}.CapMatrix

  #################################################################
  # extract the last line of the data file                         
  #################################################################
  LINE=`tail -1 ${GEOM}.CapMatrix`
  echo "${T} ${LINE}" >> ${GEOM}.W${W}_N${N}.CapVsT

done
