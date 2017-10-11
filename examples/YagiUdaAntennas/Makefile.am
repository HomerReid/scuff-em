EXTRA_DIST=                                         \
 geoFiles/Hemisphere.geo                            \
 geoFiles/PlanarYUAntenna.geo                       \
 mshFiles/Hemisphere.msh                            \
 mshFiles/PlanarYUAntenna_Medium.msh                \
 scuffgeoFiles/DipoleAntenna_Medium.scuffgeo        \
 scuffgeoFiles/YUAntenna_N2_Medium.scuffgeo         \
 scuffgeoFiles/YUAntenna_N4_Medium.scuffgeo         \
 portFiles/Dipole.ports                             \
 portcurrentFiles/Dipole.portcurrents               \
 RunScript.FieldPatterns                            \
 RunScript.SZParms

pkgdatadir = $(datadir)/scuff-em/examples/YagiUdaAntennas
nodist_pkgdata_DATA = ${EXTRA_DIST}

geodir   = $(pkgdatadir)/geoFiles
geo_DATA = geoFiles/Hemisphere.geo geoFiles/PlanarYUAntenna.geo

mshdir   = $(pkgdatadir)/mshFiles
msh_DATA = mshFiles/Hemisphere.msh mshFiles/PlanarYUAntenna_Medium.msh

portdir   = $(pkgdatadir)/portFiles
port_DATA = portFiles/Dipole.ports

scuffgeodir   = $(pkgdatadir)/scuffgeoFiles
scuffgeo_DATA = scuffgeoFiles/DipoleAntenna_Medium.scuffgeo \
                scuffgeoFiles/YUAntenna_N2_Medium.scuffgeo  \
                scuffgeoFiles/YUAntenna_N4_Medium.scuffgeo

portdir   = $(pkgdatadir)/portFiles
port_DATA = portFiles/Dipole.ports

portcurrentdir   = $(pkgdatadir)/portcurrentFiles
portcurrent_DATA = portcurrentFiles/Dipole.portcurrents
