/* -------------------------------------------------------------------
 * numpy.i
 *
 * Lisandro Dalcin <dalcinl@users.sourceforge.net>
 *
 * This SWIG file provides typemaps for NumPy
 *
 * $Id: numpy.i,v 1.1.2.2 2006/03/06 16:56:04 rodrigop Exp $
 * -----------------------------------------------------------------*/

%{
#include <numpy/arrayobject.h>
%}

%init %{
if (import_array() < 0) return;
%}


/* ---------------------------------------------------------------- */
/*  API Macros for NumPy                                            */
/* ---------------------------------------------------------------- */

/* converter functions to pass array objects from Python to C code
   as input, output and input/output */
%{
#define _PyArray_FROM_OTF(op, type, flags) \
        PyArray_FromAny((op), PyArray_DescrFromType(type), \
			0, 0, (flags), NULL)
%}


#define ARRAY_INPUT(PYOBJ, TYPENUM) \
(_PyArray_FROM_OTF((PYOBJ), (TYPENUM), CARRAY_FLAGS))

#define ARRAY_OUTPUT(PYOBJ, TYPENUM) \
(_PyArray_FROM_OTF((PYOBJ), (TYPENUM), CARRAY_FLAGS | UPDATEIFCOPY))

#define ARRAY_IO(PYOBJ, TYPENUM) \
(_PyArray_FROM_OTF((PYOBJ), (TYPENUM), CARRAY_FLAGS | UPDATEIFCOPY))


/* creation of new array objects from memory and dimension data */

%{
static PyObject *
_PyArray_NewArray(void *buffer, int type, int ndim, ...)
{
  int i;
  intp shape[MAX_DIMS];
  PyObject* array;
  va_list ap;
  va_start(ap, ndim);
  for(i=0; i<ndim; i++)
    shape[i] = va_arg(ap, intp);
  va_end(ap);
  array = PyArray_SimpleNew(ndim, shape, type);
  if (array != NULL && buffer != NULL) {
    void* data = (void*) PyArray_DATA(array);
    intp nbytes = PyArray_NBYTES(array);
    memcpy(data, buffer, nbytes);
  }
  return array;
}
%}

#define ARRAY_NEW(BUFF, NUMTYPE, ND, ...) \
(_PyArray_NewArray((void*)(BUFF), (NUMTYPE), (ND), ##__VA_ARGS__))

#define ARRAY_DEL(arr) Py_XDECREF(arr)


/* numeric type objects */

%{
static PyObject *
_PyArray_TypeObjectFromType(int type)
{
  PyArray_Descr *descr;
  PyObject *typeobj;
  descr = PyArray_DescrFromType(type);
  if (descr == NULL) return NULL;
  typeobj = (PyObject *)descr->typeobj;
  Py_XINCREF(typeobj);
  Py_DECREF(descr);
  return typeobj;
}
%}

#define ARRAY_NUMTYPEOBJ(NUMTYPE) _PyArray_TypeObjectFromType(NUMTYPE)


/* ---------------------------------------------------------------- */



/* ---------------------------------------------------------------- */
/* Error management                                                 */
/* ---------------------------------------------------------------- */

%include exception.i

%define ARRAY_arg_fail(SYMNAME, ARGNUM) 
if (SWIG_arg_fail(ARGNUM)) SWIG_fail
%enddef

%define ARRAY_exception(EXC, MESG)
SWIG_exception(EXC, "$symname(), argument number $argnum: "##MESG)
%enddef

%define ARRAY_assert(condition, message)
if (!(condition)) ARRAY_exception(SWIG_ValueError, message);
%enddef

/* ---------------------------------------------------------------- */



/* ---------------------------------------------------------------- */
/* Array checking                                                   */
/* ---------------------------------------------------------------- */

/* macros for checking array size and dimensions*/

%define ARRAY_check_size(arr, sz)
ARRAY_assert(PyArray_SIZE(arr)==(intp)(sz), "invalid array size")
%enddef

%define ARRAY_check_ndim(arr, ndim)
ARRAY_assert(PyArray_NDIM(arr)==(int)(ndim), "invalid array ndim")
%enddef

%define ARRAY_check_dim(arr, i, di)
ARRAY_assert(PyArray_DIM(arr, i)==(intp)(di), "invalid array dim")
%enddef

%define ARRAY_check_dims(array, size_t, rank, ...)
{
  size_t d;
  size_t nd         = (rank);
  size_t dims[rank] = {__VA_ARGS__};
  /* check number of dimensions */
  ARRAY_check_ndim(array, rank)
  /* check each dimension */
  for (d=0; d<nd; d++)
    ARRAY_check_dim(array, d, dims[d])
}
%enddef


/* typemaps for checking array size and dimensions */

%define ARRAY_CHECK_SIZE(TYPEMAP, SIZE)
%typemap(check) TYPEMAP 
ARRAY_check_size(array$argnum, SIZE)
%enddef

%define ARRAY_CHECK_NDIM(TYPEMAP, NDIM)
%typemap(check) TYPEMAP 
ARRAY_check_ndim(array$argnum, NDIM)
%enddef

%define ARRAY_CHECK_DIMS(TYPEMAP, SIZE_T, RANK, ...)
%typemap(check) TYPEMAP 
ARRAY_check_dims(array$argnum, SIZE_T, RANK, __VA_ARGS__)
%enddef

%define ARRAY_CHECK_DIM(TYPEMAP, I, DI)
%typemap(check) TYPEMAP 
ARRAY_check_dim(array$argnum, I, DI)
%enddef


/* ---------------------------------------------------------------- */




/* ---------------------------------------------------------------- */
/* Typecheck typemaps for overloaded C++ functions                  */
/* ---------------------------------------------------------------- */

%define ARRAY_TYPECHECK_BOOL
SWIG_TYPECHECK_BOOL_ARRAY
%enddef

%define ARRAY_TYPECHECK_INT8
SWIG_TYPECHECK_INT8_ARRAY
%enddef

%define ARRAY_TYPECHECK_INT16
SWIG_TYPECHECK_INT16_ARRAY
%enddef

%define ARRAY_TYPECHECK_INT32
SWIG_TYPECHECK_INT32_ARRAY
%enddef

%define ARRAY_TYPECHECK_INT64
SWIG_TYPECHECK_INT64_ARRAY
%enddef

%define ARRAY_TYPECHECK_INT128
SWIG_TYPECHECK_INT128_ARRAY
%enddef

%define ARRAY_TYPECHECK_FLOAT
SWIG_TYPECHECK_FLOAT_ARRAY
%enddef

%define ARRAY_TYPECHECK_DOUBLE
SWIG_TYPECHECK_DOUBLE_ARRAY
%enddef

%define ARRAY_TYPECHECK_CHAR
SWIG_TYPECHECK_CHAR_ARRAY
%enddef

%define ARRAY_TYPECHECK_STRING
SWIG_TYPECHECK_STRING_ARRAY
%enddef

%define ARRAY_TYPECHECK_OBJECT
SWIG_TYPECHECK_OBJECT_ARRAY
%enddef

%define ARRAY_TYPECHECK_SEQUENCE(TYPEMAP, PRECEDENCE)
%typemap(typecheck,precedence=PRECEDENCE) TYPEMAP
"$1 = PySequence_Check($input);";
%enddef

%define ARRAY_TYPECHECK_SEQ_OR_NUM(TYPEMAP, PRECEDENCE)
%typemap(typecheck,precedence=PRECEDENCE) TYPEMAP
"$1 = PySequence_Check($input) || PyNumber_Check($input);";
%enddef

%define ARRAY_TYPECHECK_SEQ_OR_NUM_OR_NONE(TYPEMAP, PRECEDENCE)
%typemap(typecheck, precedence=PRECEDENCE) TYPEMAP
{
  $1 = PySequence_Check($input) ||
       PyNumber_Check($input)   ||
       ($input==Py_None);
}
%enddef

/* ---------------------------------------------------------------- */



/* ---------------------------------------------------------------- */
/* ARRAY MODULE                                                     */
/* ---------------------------------------------------------------- */

%define ARRAY_ADD_MODULE(MODULE, NAME)
%init %{
{
  PyObject* module = PyImport_ImportModule(##MODULE);
  if (module == NULL) return;
  PyModule_AddObject(m, ##NAME, module);
  if (PyErr_Occurred()) return;
}
%}
%enddef

/* ---------------------------------------------------------------- */



/* ---------------------------------------------------------------- */
/* NUMERIC TYPES                                                    */
/* ---------------------------------------------------------------- */

%define ARRAY_NUMTYPE(PYNAME, MACRONAME, NUMTYPE)
#define  MACRONAME NUMTYPE /* for use in SWIG interface files */
%init %{ /* for registration in extension module */
if(PyModule_AddObject(m, #PYNAME, ARRAY_NUMTYPEOBJ(NUMTYPE))) return;
%}
%enddef

/* ---------------------------------------------------------------- */



/* ---------------------------------------------------------------- */
/* BASIC ARRAYS                                                     */
/* ---------------------------------------------------------------- */

/* raw array (this is unsafe, use with care) */

%define ARRAY_RAW(VALUE, CONVERTER, VALUE_T)
%typemap(arginit) VALUE "$1 = NULL;";
%typemap(in) VALUE (PyObject* array = NULL)
{
  if ($input != Py_None) {
    array = CONVERTER($input, VALUE_T);
    ARRAY_arg_fail($symname, $argnum);
    $1 = ($1_ltype) PyArray_DATA(array);
  }
}
%typemap(freearg) VALUE "Py_XDECREF(array$argnum);";
%enddef

%define ARRAY_RAW_CHECK_SIZE(VALUE, SIZE)
ARRAY_CHECK_SIZE(VALUE, SIZE)
%enddef

%define ARRAY_RAW_TYPECHECK(VALUE, PRECEDENCE)
ARRAY_TYPECHECK_SEQ_OR_NUM_OR_NONE(VALUE, PRECEDENCE)
%enddef


/* flat array */

%define ARRAY_FLAT(SIZE, VALUE, CONVERTER, VALUE_T)
%typemap(arginit) (SIZE, VALUE)
"$1 = 0; $2 = NULL;";
%typemap(in) (SIZE,VALUE) (PyObject* array = NULL)
{
  array = CONVERTER($input, VALUE_T);
  ARRAY_arg_fail($symname, $argnum);
  $1 = ($1_ltype) PyArray_SIZE(array);
  $2 = ($2_ltype) PyArray_DATA(array);
}
%typemap(freearg) (SIZE,VALUE)
"Py_XDECREF(array$argnum);";
%enddef

%define ARRAY_FLAT_CHECK_SIZE(SIZE, VALUE, SZ)
ARRAY_CHECK_SIZE((SIZE,VALUE), SZ)
%enddef

%define ARRAY_FLAT_TYPECHECK(SIZE, VALUE, PRECEDENCE)
ARRAY_TYPECHECK_SEQ_OR_NUM_OR_NONE((SIZE,VALUE), PRECEDENCE)
%enddef


/* ---------------------------------------------------------------- */



/* ---------------------------------------------------------------- */
/* array 1D                                                         */
/* ---------------------------------------------------------------- */

%define ARRAY_1D(SIZE, VALUE,
		 CONVERTER,VALUE_T)
%typemap(arginit) (SIZE,VALUE) 
"$1 = 0; $2 = NULL;";
%typemap(in) (SIZE,VALUE) (PyObject* array = NULL)
{
  array = CONVERTER($input, VALUE_T);
  ARRAY_arg_fail($symname,$argnum);
  $1 = ($1_ltype) PyArray_DIM(array, 0);
  $2 = ($2_ltype) PyArray_DATA(array);
}
%typemap(check) (SIZE,VALUE)
{
  if (PyArray_NDIM(array$argnum) != 1)
    ARRAY_exception(SWIG_ValueError,"array rank must be 1");
}
%typemap(freearg) (SIZE,VALUE) 
"Py_XDECREF(array$argnum);";
%enddef

/* ---------------------------------------------------------------- */


/* ---------------------------------------------------------------- */
/* array 2D                                                         */
/* ---------------------------------------------------------------- */

%define ARRAY_2D(NROW,NCOL,VALUE,
		 CONVERTER,VALUE_T)

%typemap(arginit) (NROW,NCOL,VALUE)
"$1 = 0; $2 = 0; $3 = NULL;";
%typemap(in) (NROW,NCOL,VALUE) (PyObject* array = NULL)
{
  array = CONVERTER($input, VALUE_T);
  ARRAY_arg_fail($symname,$argnum);
  $1 = ($1_ltype) PyArray_DIM(array, 0);
  $2 = ($2_ltype) PyArray_DIM(array, 1);
  $3 = ($3_ltype) PyArray_DATA(array);
}
%typemap(check) (NROW,NCOL,VALUE) 
{
  if (PyArray_NDIM(array$argnum) != 2)
    ARRAY_exception(SWIG_ValueError,"array rank must be 2");
}
%typemap(freearg) (NROW,NCOL,VALUE) 
"Py_XDECREF(array$argnum);";
%enddef

/* ---------------------------------------------------------------- */


/* ---------------------------------------------------------------- */
/* array 3D                                                         */
/* ---------------------------------------------------------------- */

%define ARRAY_3D(D1, D2, D3, VALUE,
		 CONVERTER,VALUE_T)
%typemap(arginit) (D1,D2,D3,VALUE)
"$1 = 0; $2 = 0; $3 = 0; $4 = NULL;";
%typemap(in) (D1,D2,D3,VALUE) (PyObject* array = NULL)
{
  array = CONVERTER($input, VALUE_T);
  ARRAY_arg_fail($symname, $argnum);
  $1 = ($1_ltype) PyArray_DIM(array, 0);
  $2 = ($2_ltype) PyArray_DIM(array, 1);
  $3 = ($3_ltype) PyArray_DIM(array, 2);
  $4 = ($4_ltype) PyArray_DATA(array);
}
%typemap(check) (D1,D2,D3,VALUE)
{
  if (PyArray_NDIM(array$argnum) != 3)
    ARRAY_exception(SWIG_ValueError,"array rank must be 3");
}
%typemap(freearg) (D1,D2,D3,VALUE) 
"Py_XDECREF(array$argnum);";
%enddef

/* ---------------------------------------------------------------- */


/* ---------------------------------------------------------------- */
/* array ND                                                         */
/* ---------------------------------------------------------------- */

%define ARRAY_ND(RANK, DIMS, VALUE,
		 CONVERTER,VALUE_T)
%typemap(arginit) (RANK,DIMS,VALUE)
"$1 = 0; $2 = NULL; $3 = NULL;";
%typemap(in) (RANK,DIMS,VALUE) (PyObject* array = NULL)
{
  /* get the correct Numarray (C compatible) */
  array = CONVERTER($input, VALUE_T);
  ARRAY_arg_fail($symname, $argnum);
  /* get Numarray's dim & data */
  $1 = ($1_ltype) PyArray_NDIM(array);
  $2 = ($2_ltype) PyMem_New($*2_ltype, $1);
  if(!$2) ARRAY_exception(SWIG_MemoryError, "");
  {$1_ltype d; for (d=0; d<$1; d++)
    $2[d] = ($*2_ltype) PyArray_DIM(array, d);}
  $3 = ($3_ltype) PyArray_DATA(array);
}
%typemap(freearg) (RANK,DIMS,VALUE) 
"PyMem_Del($2); Py_XDECREF(array$argnum);";
%enddef

/* ---------------------------------------------------------------- */




/* ---------------------------------------------------------------- */
/* Array pair from a sequence                                       */
/* ---------------------------------------------------------------- */

%define ARRAY_PAIR(SIZE, VALUE1, VALUE2,
		   CONVERTER1, VALUE1_T,
		   CONVERTER2, VALUE2_T)
%typemap(arginit) (SIZE, VALUE1, VALUE2)
"$1 = 0; $2 = NULL; $3 = NULL; "
%typemap(in) (SIZE,VALUE1,VALUE2)
     (PyObject* oarray1 = NULL, PyObject* array1 = NULL,
      PyObject* oarray2 = NULL, PyObject* array2 = NULL)
{
  /* parse input sequence */
  PyObject* targs = PySequence_Tuple($input);
  ARRAY_arg_fail($symname, $argnum);
  PyArg_UnpackTuple(targs, "", 2, 2,
		    &oarray1, &oarray2);
  Py_DECREF(targs);
  ARRAY_arg_fail($symname, $argnum);
  /* obtain array object */
  array1 = CONVERTER1(oarray1, VALUE1_T);
  ARRAY_arg_fail($symname, $argnum);
  array2 = CONVERTER2(oarray2, VALUE2_T);
  ARRAY_arg_fail($symname, $argnum);
  /* extract array data */
  $1 = ($1_ltype) PyArray_SIZE(array1);
  $2 = ($2_ltype) PyArray_DATA(array1);
  $3 = ($3_ltype) PyArray_DATA(array2);
}
%typemap(freearg) (SIZE,VALUE1,VALUE2)
"Py_XDECREF(array1$argnum); Py_XDECREF(array2$argnum);";
%enddef

%define ARRAY_PAIR_CHECK_SIZE(SIZE, VALUE1, VALUE2)
%typemap(check) (SIZE,VALUE1,VALUE2)
{ ARRAY_check_size(array2$argnum, $1) }
%enddef

%define ARRAY_PAIR_TYPECHECK(SIZE, VALUE1, VALUE2, PRECEDENCE)
%typemap(typecheck, precedence=PRECEDENCE) (SIZE,VALUE1,VALUE2)
{
  $1 = (PySequence_Check($input) &&
	PySequence_Size($input)==2 )? 1 : 0;
}
%enddef


/* ---------------------------------------------------------------- */


/* ---------------------------------------------------------------- */
/* Array triad from a sequence                                      */
/* ---------------------------------------------------------------- */

%define ARRAY_TRIAD(SIZE, VALUE1, VALUE2, VALUE3,
		    CONVERTER1, VALUE1_T,
		    CONVERTER2, VALUE2_T,
		    CONVERTER3, VALUE3_T)
%typemap(arginit) (SIZE,VALUE1,VALUE2,VALUE3)
"$1 = 0; $2 = NULL; $3 = NULL; $4 = NULL;"
%typemap(in) (SIZE,VALUE1,VALUE2,VALUE3)
     (PyObject* oarray1 = NULL, PyObject* array1 = NULL,
      PyObject* oarray2 = NULL, PyObject* array2 = NULL,
      PyObject* oarray3 = NULL, PyObject* array3 = NULL)
{
  /* parse input sequence */
  PyObject* targs = PySequence_Tuple($input);
  ARRAY_arg_fail($symname, $argnum);
  PyArg_UnpackTuple(targs, "", 3, 3,
		    &oarray1, &oarray2, &oarray3);
  Py_DECREF(targs);
  ARRAY_arg_fail($symname, $argnum);
  /* obtain array object */
  array1 = CONVERTER1(oarray1, VALUE1_T);
  ARRAY_arg_fail($symname, $argnum);
  array2 = CONVERTER2(oarray2, VALUE2_T);
  ARRAY_arg_fail($symname, $argnum);
  array3 = CONVERTER3(oarray3, VALUE3_T);
  ARRAY_arg_fail($symname, $argnum);
  /* extract array data */
  $1 = ($1_ltype) PyArray_SIZE(array1);
  $2 = ($2_ltype) PyArray_DATA(array1);
  $3 = ($3_ltype) PyArray_DATA(array2);
  $4 = ($4_ltype) PyArray_DATA(array3);
}
%typemap(freearg) (SIZE,VALUE1,VALUE2,VALUE3)
{
  Py_XDECREF(array1$argnum);
  Py_XDECREF(array2$argnum);
  Py_XDECREF(array3$argnum);
}
%enddef

%define ARRAY_TRIAD_CHECK_SIZE(SIZE,VALUE1,VALUE2,VALUE3)
%typemap(check) (SIZE,VALUE1,VALUE2,VALUE3)
{ 
  ARRAY_check_size(array2$argnum, $1)
  ARRAY_check_size(array3$argnum, $1)
}
%enddef

%define ARRAY_TRIAD_TYPECHECK(SIZE, VALUE1, VALUE2, VALUE3, PRECEDENCE)
%typemap(typecheck, precedence=PRECEDENCE) (SIZE,VALUE1,VALUE2,VALUE3)
{
  $1 = (PySequence_Check($input) &&
	PySequence_Size($input)==3 )? 1 : 0;
}
%enddef


/* ---------------------------------------------------------------- */



/* ---------------------------------------------------------------- */
/* FEM element-level arrays                                         */
/* ---------------------------------------------------------------- */

/* elemental vector */

%define ARRAY_ELEM_VEC(SIZE, INDEX, VALUE,
		       IDX_CONVERTER, IDX_T,
		       IDX_CONVERTER, VAL_T)
ARRAY_PAIR(SIZE, INDEX, VALUE,
	   IDX_CONVERTER, IDX_T,
	   IDX_CONVERTER, VAL_T)
ARRAY_PAIR_CHECK_SAME_SIZE(SIZE, INDEX, VALUE)
%enddef

%define ARRAY_ELEM_VEC_TYPECHECK(SIZE, INDEX, VALUE, PRECEDENCE)
ARRAY_PAIR_TYPECHECK(SIZE, INDEX, VALUE, PRECEDENCE)
%enddef


/* elemental matrix */

%define ARRAY_ELEM_MAT(NROW,ROW,NCOL,COL,VALUE,
		       ROW_CONVERTER, ROW_T,
		       COL_CONVERTER, COL_T,
		       VAL_CONVERTER, VAL_T)
%typemap(arginit) (NROW,ROW,NCOL,COL,VALUE)
"$1 = 0; $2 = NULL; $3 = 0; $4 = NULL; $5 = NULL; "
%typemap(in) (NROW,ROW,NCOL,COL,VALUE)
     (PyObject* oarray1 = NULL, PyObject* array1 = NULL,
      PyObject* oarray2 = NULL, PyObject* array2 = NULL,
      PyObject* oarray3 = NULL, PyObject* array3 = NULL)
{
  /* parse input sequence */
  PyObject* targs = PySequence_Tuple($input);
  ARRAY_arg_fail($symname, $argnum);
  PyArg_UnpackTuple(targs, "", 3, 3,
		    &oarray1, &oarray2, &oarray3);
  Py_DECREF(targs);
  ARRAY_arg_fail($symname, $argnum);
  /* obtain array objects */
  array1 = ROW_CONVERTER(oarray1, ROW_T);
  ARRAY_arg_fail($symname, $argnum);
  array2 = COL_CONVERTER(oarray2, COL_T);
  ARRAY_arg_fail($symname, $argnum);
  array3 = VAL_CONVERTER(oarray3, VAL_T);
  ARRAY_arg_fail($symname, $argnum);
  /* extract data from arrays */
  $1 = ($1_ltype) PyArray_SIZE(array1);
  $2 = ($2_ltype) PyArray_DATA(array1);
  $3 = ($3_ltype) PyArray_SIZE(array2);
  $4 = ($4_ltype) PyArray_DATA(array2);
  $5 = ($5_ltype) PyArray_DATA(array3);
}
%typemap(check) (NROW,ROW,NCOL,COL,VALUE)
{ 
  ARRAY_check_size(array3$argnum, $1*$3)
}
%typemap(freearg) (NROW,ROW,NCOL,COL,VALUE)
{
  Py_XDECREF(array1$argnum);
  Py_XDECREF(array2$argnum);
  Py_XDECREF(array3$argnum);
}
%enddef

%define ARRAY_ELEM_MAT_TYPECHECK(NROW,ROW,NCOL,COL,VALUE, PRECEDENCE)
%typemap(typecheck, precedence=PRECEDENCE) (NROW,ROW,NCOL,COL,VALUE)
{
  $1 = (PySequence_Check($input) && 
	PySequence_Size($input)==3 )? 1 : 0;
}
%enddef

/* ---------------------------------------------------------------- */



/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* OUTPUT TYPEMAPS                                                  */
/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */


/* ---------------------------------------------------------------- */
/* array 1D                                                         */
/* ---------------------------------------------------------------- */

%define ARRAY_1D_NEW(SIZE, VALUE, VALUE_T)
%typemap(in,numinputs=0) (SIZE,VALUE)
     ($*1_ltype size, 
      $*2_ltype value)
{ 
  size  = ($*1_ltype) 0; $1 = &size; 
  value = ($*2_ltype) 0; $2 = &value;
}
%typemap(argout) (SIZE, VALUE)
{
  PyObject* o = NULL;
  o = ARRAY_NEW(*$2, VALUE_T, 1, *$1);
  ARRAY_arg_fail($symname, $argnum);
  %append_output(o);
}
%enddef

%define ARRAY_1D_FREEARG(SIZE, VALUE, FREE)
%typemap(freearg) (SIZE,VALUE)
"if (*$2) FREE(*$2);";
%enddef

/* ---------------------------------------------------------------- */


/* ---------------------------------------------------------------- */
/* array 2D                                                         */
/* ---------------------------------------------------------------- */

%define ARRAY_2D_NEW(ROW, COL, VALUE, VALUE_T)

%typemap(in,numinputs=0) (ROW,COL,VALUE)
     ($*1_ltype row,
      $*2_ltype col,
      $*3_ltype val)
{ 
  row = ($*1_ltype) 0; $1 = &row;
  col = ($*2_ltype) 0; $2 = &col;
  val = ($*3_ltype) 0; $3 = &val;
}

%typemap(argout) (ROW, COL, VALUE)
{
  PyObject* o = NULL;
  o = ARRAY_NEW(*$3, VALUE_T, 2, *$1, *$2);
  ARRAY_arg_fail($symname, $argnum);
  %append_output(o);
}

%enddef

%define ARRAY_2D_FREEARG(ROW, COL, VALUE, FREE)
%typemap(freearg) (ROW,COL,VALUE)
"if (*$3) FREE(*$3);";
%enddef


/* ---------------------------------------------------------------- */


/* ---------------------------------------------------------------- */
/* array 3D                                                         */
/* ---------------------------------------------------------------- */

%define ARRAY_3D_NEW(D1, D2, D3, VALUE, VALUE_T)

%typemap(in, numinputs=0) (D1, D2, D3, VALUE)
     ($*1_ltype d1,
      $*2_ltype d2,
      $*3_ltype d3,
      $*4_ltype vl)
{ 
  d1 = ($*1_ltype) 0; $1 = &d1; 
  d2 = ($*2_ltype) 0; $2 = &d2; 
  d3 = ($*3_ltype) 0; $3 = &d3; 
  vl = ($*4_ltype) 0; $4 = &vl;
}

%typemap(argout) (D1, D2, D3, VALUE)
{
  PyObject* o = NULL;
  o = ARRAY_NEW(*$4, VALUE_T, 3, *$1, *$2,*$3);
  ARRAY_arg_fail($symname, $argnum);
  %append_output(o);
}

%enddef

%define ARRAY_3D_FREEARG(D1, D2, D3, VALUE, FREE)
%typemap(freearg) (D1,D2,D3,VALUE)
"if (*$4) FREE(*$4);";
%enddef

/* ---------------------------------------------------------------- */





/*
 * Local Variables:
 * mode: C
 * End:
 */
