// -*- c++ -*-
// $Id: forward.h,v 1.1.2.6 2006/03/30 15:40:05 rodrigop Exp $

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
class Dofmap;

typedef TextHashTable OptionTable;
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
class DofMap;

typedef ::OptionTable OptionTable;

PYPF_NAMESPACE_END

#endif // PYPF_FORWARD_H
