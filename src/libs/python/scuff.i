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

#ifdef SWIGPYTHON
%include "scuff-python.i"
#endif

%include "libscuff.h"
