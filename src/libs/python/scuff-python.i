// -*- C++ -*-

%{
#define SWIG_FILE_WITH_INIT
#define array_stride(a,i)        (((PyArrayObject *)a)->strides[i])
%}
%include "numpy.i"
%init %{
  import_array();
%}
%fragment("NumPy_Fragments");
%fragment("NumPy_Macros");
%numpy_typemaps(double, NPY_DOUBLE, size_t)
%numpy_typemaps(cdouble, NPY_CDOUBLE, size_t)

//////////////////////////////////////////////////////////////////////////////
// Wrapper for EHFuncType

// wrap Python function EH = f(R)
%{
static void EHFunc_python(double *R, void *f, cdouble *EH) {
  npy_intp sz3 = 3, stride1 = sizeof(double);
  PyObject *Rpy = PyArray_New(&PyArray_Type, 1, &sz3, NPY_DOUBLE, &stride1,
                              R, // not NPY_WRITEABLE, effectively const
                              0, NPY_C_CONTIGUOUS | NPY_ALIGNED, NULL);
  PyObject *arglist = Py_BuildValue("O", Rpy);
  PyObject *result = PyEval_CallObject((PyObject *) f, arglist);
  PyArrayObject *EHpy;
  int is_new;
  Py_DECREF(arglist);
  Py_DECREF(Rpy);
  if (PyErr_Occurred()) {
    ErrExit("python exception handling not supported in libscuff");
  }
  else if ((EHpy = obj_to_array_contiguous_allow_conversion(result,
                                                            NPY_CDOUBLE,
                                                            &is_new))) {
    if (array_numdims(EHpy) != 1 || array_size(EHpy,0) != 0)
      ErrExit("6-component EH field must be returned to libscuff");
    memcpy(EH, array_data(EHpy), sizeof(cdouble) * 6);
    if (is_new) Py_DECREF((PyObject *) EHpy);
  }
  else {
    ErrExit("invalid EH field returned to libscuff");
  }
  Py_XDECREF(result);
}
%}

%typemap(in)(scuff::EHFuncType EHFunc, void *EHFuncUD) {
  $1 = EHFunc_python;
  Py_INCREF($input);
  $2 = (void *) $input;
}
%typemap(freearg)(scuff::EHFuncType EHFunc, void *EHFuncUD) {
  Py_XDECREF((PyObject *) $2);
}
%typecheck(SWIG_TYPECHECK_POINTER)(scuff::EHFuncType EHFunc, void *EHFuncUD) {
  $1 = PyCallable_Check($input);
}

//////////////////////////////////////////////////////////////////////////////
// HMatrix and HVector typemaps: these are mapped directly onto
// NumPy matrices, ideally with no data copying.
//
// TODO: support opaque pointers for packed (symmetric) HMatrices.

%typemap(out)(HVector *) {
  npy_intp sz = npy_intp($1->N);
  $result = (PyObject*) PyArray_SimpleNewFromData
    (1, &sz, $1->RealComplex ? NPY_CDOUBLE : NPY_DOUBLE,
     $1->RealComplex ? (void*) $1->ZV : (void*) $1->DV);
  $1->ZV = 0; $1->DV = 0; // prevent deallocation
  delete $1;
}

%{
  static bool is_HVector(PyObject *o) {
    return is_array(o) && array_numdims(o) == 1
      && array_is_contiguous(o) && array_is_native(o)
      && PyArray_ISALIGNED(o)
      && (array_type(o) == NPY_DOUBLE || array_type(o) == NPY_CDOUBLE);
  }
%}

%typemap(in)(HVector *) {
  if (is_HVector($input)) {
    // hack to create HVector without making a copy of the data
    $1 = new HVector(0, array_type($input) == NPY_DOUBLE 
		        ? LHM_REAL : LHM_COMPLEX);
    if (array_type($input) == NPY_DOUBLE)
      $1->DV = (double*) array_data($input);
    else
      $1->ZV = (cdouble*) array_data($input);
    $1->N = array_size($input, 0);
  }
  else {
    $1 = NULL; // TODO: allow conversions (with a copy) from other types?
    ErrExit("invalid HVector argument - expecting numpy.array");
  }
}
%typemap(freearg)(HVector *) {
  $1->DV = 0; $1->ZV = 0; delete $1; // hack to avoid deallocating Python data
}
%typecheck(SWIG_TYPECHECK_POINTER)(HVector *) {
  $1 = is_HVector($input);
}

%typemap(out)(HMatrix *) {
  if ($1->StorageType != LHM_NORMAL) { // TODO: return opaque HMatrix* ?
    HMatrix *copy_$1 = new HMatrix($1->NR, $1->NC, $1->RealComplex);
    copy_$1->Copy($1);
    delete $1;
    $1 = copy_$1;
  }
  npy_intp dims[2];
  dims[0] = npy_intp($1->NR);
  dims[1] = npy_intp($1->NR);
  $result = (PyObject*) PyArray_New // TODO: how to return numpy.matrix?
    (&PyArray_Type, 2, dims, $1->RealComplex ? NPY_CDOUBLE : NPY_DOUBLE,
     NULL, $1->RealComplex ? (void*) $1->ZM : (void*) $1->DM, 0,
     NPY_FARRAY, NULL);
  $1->ZM = 0; $1->DM = 0; // prevent deallocation
  delete $1;
}

%{
  static bool is_HMatrix(PyObject *o) {
    return is_array(o) && array_numdims(o) == 2
      && array_is_contiguous(o) && array_is_native(o)
      && PyArray_ISALIGNED(o)
      && (array_type(o) == NPY_DOUBLE || array_type(o) == NPY_CDOUBLE);
  }
%}

%typemap(in)(HMatrix *) {
  if (is_HMatrix($input)) {
    // hack to create HMatrix without making a copy of the data
    $1 = new HMatrix(0, 0, array_type($input) == NPY_DOUBLE 
		        ? LHM_REAL : LHM_COMPLEX);
    if (array_type($input) == NPY_DOUBLE)
      $1->DM = (double*) array_data($input);
    else
      $1->ZM = (cdouble*) array_data($input);
    if (array_is_fortran($input)) {
      $1->NR = array_size($input, 0);
      $1->NC = array_size($input, 1);
    }
    else { // C-ordered: transposed dimensions (FIXME: transpose the data?)
      $1->NR = array_size($input, 1);
      $1->NC = array_size($input, 0);
    }
  }
  else {
    $1 = NULL; // TODO: allow conversions (with a copy) from other types?
    ErrExit("invalid HMatrix argument - expecting numpy.array");
  }
}
%typemap(freearg)(HMatrix *) {
  $1->DM = 0; $1->ZM = 0; delete $1; // hack to avoid deallocating Python data
}
%typecheck(SWIG_TYPECHECK_POINTER)(HMatrix *) {
  $1 = is_HMatrix($input);
}
