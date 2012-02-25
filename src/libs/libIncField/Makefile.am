lib_LIBRARIES = libIncField.a

include_HEADERS = libIncField.h

libIncField_a_SOURCES = 	\
 IncFieldData.o			\
 PlaneWave.o    		\
 PointSource.o	 		\
 libIncField.h

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libhrutil    \
              -I$(top_srcdir)/src/libs/libMatProp   \
              -I$(top_srcdir)/src/libs/libhmat      \
              -I$(top_srcdir)/src/libs/libSGJC      \
              -I$(top_srcdir)/src/libs/libIncField  \
              -I$(top_srcdir)/src/libs/libMDInterp  \
              -I$(top_srcdir)/src/libs/libTriInt
