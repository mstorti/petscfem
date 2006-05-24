// $Id: forward.h,v 1.1.2.8 2006/05/24 20:57:53 dalcinl Exp $

#ifndef PYPF_FORWARD_H
#define PYPF_FORWARD_H


// PETSc
// -----

typedef struct _p_Vec* Vec;
typedef struct _p_Mat* Mat;


// PETScFEM
// --------

class TextHashTable;
class Nodedata;
class Elemset;
class Mesh;
class Constraint;
class Amplitude;
class Dofmap;

typedef TextHashTable Options;
typedef Nodedata      Nodeset;
typedef Dofmap        DofMap;


// PETScFEM for Python
// -------------------

#include "namespace.h"

PYPF_NAMESPACE_BEGIN

class Options;
class Object;
class Nodeset;
class Elemset;
class Mesh;
class Amplitude;
class DofMap;
class Domain;

class Problem;

class Application;

PYPF_NAMESPACE_END

#endif // PYPF_FORWARD_H

// Local Variables:
// mode: C++
// End:
