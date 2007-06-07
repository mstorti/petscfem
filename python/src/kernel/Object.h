// $Id$

#ifndef PF4PY_OBJECT_H
#define PF4PY_OBJECT_H

#include "petscfem4py.h"

#include "RefCnt.h"

PF4PY_NAMESPACE_BEGIN

class Object
#if !defined(SWIG)
  : public RefCnt
#endif
{

public:
  virtual ~Object() = 0;
  Object() { };
  Object(const Object&) { };

public:
  bool operator==(const Object& obj) const { return this == &obj; }
  bool operator!=(const Object& obj) const { return this != &obj; }

};

PF4PY_NAMESPACE_END

#endif // PF4PY_OBJECT_H

// Local Variables:
// mode: C++
// End:
