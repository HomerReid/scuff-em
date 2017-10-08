EXTRA_DIST=				        \
 geoFiles/Hemisphere.geo			\
 geoFiles/PlanarYUAntenna.geo  			\
 mshFiles/Hemisphere.msh       			\
 mshFiles/PlanarYUAntenna_Medium.msh 		\
 portFiles/Dipole.ports        			\
 portcurrentFiles/Dipole.portCurrents 		\
 scuffgeoFiles/DipoleAntenna_Medium.scuffgeo 	\
 scuffgeoFiles/YUAntenna_N2_Medium.scuffgeo 	\
 scuffgeoFiles/YUAntenna_N4_Medium.scuffgeo 	\
 RunScript.FieldPatterns                    	\
 RunScript.SZParms

pkgdatadir = $(datadir)/scuff-em/examples/YagiUdaAntennas
nodist_pkgdata_DATA = ${EXTRA_DIST}

geodir   = $(pkgdatadir)/geoFiles
geo_DATA = geoFiles/Hemisphere.geo geoFiles/PlanarYUAntenna.geo

mshdir   = $(pkgdatadir)/mshFiles
msh_DATA = mshFiles/Hemisphere.msh mshFiles/PlanarYUAntenna_Medium.msh

portdir   = $(pkgdatadir)/portFiles
port_DATA = portFiles/Dipole.ports

portcurrentdir   = $(pkgdatadir)/portcurrentFiles
portcurrent_DATA = portcurrentFiles/Dipole.portCurrents

scuffgeodir   = $(pkgdatadir)/scuffgeoFiles
scuffgeo_DATA = scuffgeoFiles/DipoleAntenna_Medium.scuffgeo \
                scuffgeoFiles/YUAntenna_N2_Medium.scuffgeo  \
                scuffgeoFiles/YUAntenna_N4_Medium.scuffgeo
