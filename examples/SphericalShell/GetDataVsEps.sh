#!/bin/bash

# This script runs scuff-scatter to compute the E-field at 
# the origin for various different values of the shell 
# permittivity.

cat EpsValues | while read EPS
do

 # 
 # create a copy of the .scuffgeo file with EPS_10 replaced
 # by the current value of the permittivity
 # 
 sed "s/EPS_10/EPS_${EPS}/" SphericalShell.scuffgeo > temp.scuffgeo

 # 
 # run scuff-scatter to get the E-field at the origin for this 
 # value of epsilon
 # 
 /bin/rm -f CenterPoint.total  
 /bin/rm -f CenterPoint.scattered
 scuff-scatter --geometry temp.scuffgeo --EPFile CenterPoint < Args 

 # extract the z-component of the field from the output file
 EZ=`awk '{print $8}' CenterPoint.total`
 echo "${EPS} ${EZ}" >> EzVsEps.dat

done
