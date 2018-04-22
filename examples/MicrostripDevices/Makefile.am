EXTRA_DIST=				\
 geoFiles/Rectangle.geo			\
 geoFiles/EdgeFedAntenna.geo		\
 geoFiles/PointFedAntenna.geo		\
 mshFiles/Center_436.msh		\
 mshFiles/Center_237.msh		\
 mshFiles/Center_98.msh			\
 mshFiles/Flank_912.msh              	\
 mshFiles/Flank_504.msh              	\
 mshFiles/Flank_216.msh              	\
 mshFiles/EFAntenna_L0_Coarse.msh       \
 mshFiles/EFAntenna_L0_Fine.msh         \
 mshFiles/EFAntenna_L8_Coarse.msh       \
 mshFiles/EFAntenna_L8_Fine.msh         \
 mshFiles/PointFedAntenna_Coarse.msh    \
 mshFiles/PointFedAntenna_Fine.msh      \
 scuffgeoFiles/CPW_Coarse.scuffgeo   	\
 scuffgeoFiles/CPW_Medium.scuffgeo   	\
 scuffgeoFiles/CPW_Fine.scuffgeo     	\
 scuffgeoFiles/CPW_Fine.scuffgeo     	\
 portFiles/CPW.ports			\
 portFiles/EFAntenna.ports		\
 portFiles/CoupledAntennas.ports	\
 transFiles/HPlane.trans        	\
 transFiles/EPlane.trans

pkgdatadir = $(datadir)/scuff-em/examples/MicrostripDevices
nodist_pkgdata_DATA = ${EXTRA_DIST}

geodir   = $(pkgdatadir)/geoFiles
geo_DATA = geoFiles/Rectangle.geo       		\
           geoFiles/EdgeFedAntenna.geo			\
           geoFiles/PointFedAntenna.geo

mshdir   = $(pkgdatadir)/mshFiles
msh_DATA = mshFiles/Center_436.msh 			\
           mshFiles/Center_237.msh 			\
           mshFiles/Center_98.msh  			\
           mshFiles/Flank_912.msh			\
           mshFiles/Flank_504.msh                 	\
           mshFiles/Flank_216.msh                 	\
           mshFiles/EFAntenna_L0_Coarse.msh       	\
           mshFiles/EFAntenna_L0_Fine.msh         	\
           mshFiles/EFAntenna_L8_Coarse.msh       	\
           mshFiles/EFAntenna_L8_Fine.msh		\
           mshFiles/PointFedAntenna_Coarse.msh		\
           mshFiles/PointFedAntenna_Fine.msh

scuffgeodir   = $(pkgdatadir)/scuffgeoFiles
scuffgeo_DATA = scuffgeoFiles/CPW_Coarse.scuffgeo		\
                scuffgeoFiles/CPW_Medium.scuffgeo		\
                scuffgeoFiles/CPW_Fine.scuffgeo			\
                scuffgeoFiles/EFAntenna_L0_Coarse.scuffgeo	\
                scuffgeoFiles/EFAntenna_L0_Fine.scuffgeo	\
                scuffgeoFiles/EFAntenna_L8_Coarse.scuffgeo	\
                scuffgeoFiles/EFAntenna_L8_Fine.scuffgeo	\
                scuffgeoFiles/CoupledAntennas_Coarse.scuffgeo	\
                scuffgeoFiles/CoupledAntennas_Fine.scuffgeo

portdir   = $(pkgdatadir)/portFiles
port_DATA = portFiles/CPW.ports 				\
            portFiles/EFAntenna.ports				\
            portFiles/CoupledAntennas.ports 

transdir   = $(pkgdatadir)/transFiles
trans_DATA = transFiles/EPlane.trans				\
             transFiles/HPlane.trans
