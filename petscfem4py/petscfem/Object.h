// $Id: Object.h,v 1.1.2.9 2006/08/22 22:10:43 dalcinl Exp $

#ifndef PYPF_OBJECT_H
#define PYPF_OBJECT_H

#include <string>
#include <map>
#include <mpi.h>
#include "petscfem4py.h"
#include "Options.h"

PYPF_NAMESPACE_BEGIN

class Object
#if !defined(SWIG)
  : public RefCounter
#endif
{

protected:
  Object(MPI_Comm comm);
public:
  virtual ~Object() = 0;
  Object();
  Object(const Object&);

  // object comparisons
public:
  bool operator==(const Object& obj) const { return this == &obj; }
  bool operator!=(const Object& obj) const { return this != &obj; }

  // MPI communicator
protected:
  MPI_Comm comm;
public:
  MPI_Comm getComm() const;
  void     setComm(MPI_Comm comm);

  // options management
public:
  Options options;

};

PYPF_NAMESPACE_END

#endif // PYPF_OBJECT_H

// Local Variables:
// mode: C++
// End:
