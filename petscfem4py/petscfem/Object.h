// $Id: Object.h,v 1.1.2.6 2006/04/27 19:09:17 rodrigop Exp $

#ifndef PYPF_OBJECT_H
#define PYPF_OBJECT_H

#include <string>
#include <map>
#include <mpi.h>
#include "petscfem4py.h"
#include "Options.h"
#include "Comm.h"

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

  // reference count
public:
  int getRefCount() const { return getref(); }

  // object comparisons
public:
  bool operator==(const Object& obj) const { return this == &obj; }
  bool operator!=(const Object& obj) const { return this != &obj; }

  // communicator
protected:
  MPI_Comm comm;
public:
  Comm getComm() const;
  void setComm(Comm& comm);
  void setComm(MPI_Comm comm);

  // options management
protected:
  Options options;
public:
  bool        hasOption (const std::string& key) const;
  std::string getOption (const std::string& key) const;
  void        setOption (const std::string& key, const std::string& value);
  std::map<std::string,std::string> getOptions() const;
  void setOptions(const std::map<std::string,std::string>& options);
  void addOptions(const std::map<std::string,std::string>& options);
  void delOptions();

};


PYPF_NAMESPACE_END

#endif // PYPF_OBJECT_H

// Local Variables:
// mode: C++
// End:
