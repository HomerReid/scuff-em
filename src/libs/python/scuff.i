// -*- C++ -*-

%define DOCSTRING
"SCUFF-EM is a free open-source software suite for boundary-element
analysis of problems in computational electromagnetism and related
fields."
%enddef

%module(docstring=DOCSTRING) scuff
%{
#include "libscuff.h"
using namespace scuff;
%}

//////////////////////////////////////////////////////////////////////////////
// Need to explicitly tell SWIG which functions allocate new objects
// that the caller is responsible for freeing.
// FIXME: SWIG seems to be ignoring these directives(?)

%newobject HVector::Copy();
%newobject HMatrix::Copy();
%newobject CopyHMatrixUnpacked;
%newobject scuff::RWGGeometry::AllocateBEMMatrix;
%newobject scuff::RWGGeometry::AllocateDMDVMatrix;
%newobject scuff::RWGGeometry::AllocateRHSVector;
%typemap(newfree) GTransformation * "free($1);"
%newobject scuff::CreateGTransformation;

//////////////////////////////////////////////////////////////////////////////

#ifdef SWIGPYTHON
%include "scuff-python.i"
#endif

//////////////////////////////////////////////////////////////////////////////

%include "libhmat.h"
%include "libMatProp.h"
%include "libIncField.h"
%include "GTransformation.h"
%include "libscuff.h"
