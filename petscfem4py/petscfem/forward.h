// $Id: forward.h,v 1.1.2.7 2006/04/27 19:09:17 rodrigop Exp $

#ifndef PYPF_FORWARD_H
#define PYPF_FORWARD_H


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

PYPF_NAMESPACE_END

#endif // PYPF_FORWARD_H

// Local Variables:
// mode: C++
// End:
