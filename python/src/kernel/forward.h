// $Id$

#ifndef PF4PY_FORWARD_H
#define PF4PY_FORWARD_H

// PETSc
// -----

typedef struct _p_Vec*        Vec;
typedef struct _p_VecScatter* VecScatter;
typedef struct _p_Mat*        Mat;


// PETScFEM
// --------

class TextHashTable;
class Nodedata;
class Elemset;
class Mesh;
class Constraint;
class Amplitude;
class Dofmap;
class arg_list;
class State;
class Time;


// PETScFEM for Python
// -------------------

#include "namespace.h"

PF4PY_NAMESPACE_BEGIN

class Comm;
class Object;
class Options;
class Elemset;
class Mesh;
class Amplitude;
class Dofset;
class AppCtx;
class Domain;

PF4PY_NAMESPACE_END

#endif // PF4PY_FORWARD_H

// Local Variables:
// mode: C++
// End:
