// -*- C++ -*-

%define DOCSTRING
"SCUFF-EM is a free open-source software suite for boundary-element
analysis of problems in computational electromagnetism and related
fields."
%enddef

%module(docstring=DOCSTRING) scuff

%include "scuff-python.i"

%{
#include "libscuff.h"
#include "scuffSolver.h"
#include "StaticSolver.h"
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
%newobject scuff::StaticSolver::AllocateBEMMatrix;
%newobject scuff::StaticSolver::AllocateRHSVector;
%typemap(newfree) GTransformation * "free($1);"
%newobject scuff::CreateGTransformation;

// SWIG has difficulty with functions that take va_list arguments
%ignore vsnprintfEC;

// these functions are defined in two different libraries
// but we do not need them in python anyway
%ignore VecScale;
%ignore VecPlusEquals;

//////////////////////////////////////////////////////////////////////////////

%include "libhrutil.h"
%include "libhmat.h"
%include "libMatProp.h"
%include "libIncField.h"
%include "GTransformation.h"
%include "libscuff.h"
%include "StaticSolver.h"
%include "scuffSolver.h"
