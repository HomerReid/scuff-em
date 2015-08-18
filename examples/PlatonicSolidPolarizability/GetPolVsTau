#!/bin/bash

##################################################
# specify fixed file names, paths, etc
##################################################
CODE=scuff-static
GEOTEMPLATE=Template.scuffgeo
GEOFILE=Temp.scuffgeo
POLFILE=Temp.pol
EPSFILE=EpsFile

##################################################
# loop over geometry meshes
##################################################
SHAPELIST="Sphere_498 Tetrahedron_615 Octahedron_456 Sphere_924 Octahedron_975 Tetrahedron_1179"
for SHAPE in ${SHAPELIST}
do

  echo "Computing polarizability vs. Tau for ${SHAPE}..."
  OUTFILE=${SHAPE}.AlphaVsTau
  /bin/rm -rf ${OUTFILE}

  ##################################################
  # loop over (EpsOut, EpsIn) pairs
  ##################################################
  cat EPSFILE | while read EPSOUT EPSIN
  do

    echo " Eps{Out,In} = {${EPSOUT}, ${EPSIN}}... "
    ##################################################
    # create a .scuffgeo file for this (Shape,EpsOut,EpsIn)
    ##################################################
    /bin/cp ${GEOTEMPLATE} ${GEOFILE}
    sed -i "s/EPSIN/${EPSOUT}/"  ${GEOFILE}
    sed -i "s/EPSOUT/${EPSIN}/"  ${GEOFILE}
    sed -i "s/SHAPE/${SHAPE}/"   ${GEOFILE}
    
    ##################################################
    # run scuff-static to get the polarizability
    ##################################################
    /bin/rm -rf ${DATAFILE}
    ${CODE} --geometry ${GEOFILE} --polfile ${POLFILE} < /dev/null

    ##################################################
    # grab the polarizability from the output file written by 
    # scuff-static and dump it together with (EpsOut,EpsIn)
    # data to the overall output file for this shape
    ##################################################
    printf "${EPSOUT} ${EPSIN} " >> ${OUTFILE}
    tail -1 ${POLFILE} >> ${OUTFILE}
    echo >> ${OUTFILE}

  done

  echo "Polarizability vs Tau written to ${OUTFILE}."
done

echo "Thank you for your support."
