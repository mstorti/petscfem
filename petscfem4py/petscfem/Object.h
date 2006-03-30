// -*- c++ -*-
// $Id: Object.h,v 1.1.2.5 2006/03/30 15:40:05 rodrigop Exp $

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
  MPI_Comm& getComm() const;

  // options management
protected:
  Options options;
public:
  typedef std::string             string;
  typedef std::map<string,string> Table;
  bool   hasOption(const string& key) const;
  string getOption(const string& key) const;
  void   setOption(const string& key, const string& value);
  Table getOptions() const;
  void  setOptions(const Table&);
  void  addOptions(const Table&);

};


PYPF_NAMESPACE_END

#endif // PYPF_OBJECT_H
