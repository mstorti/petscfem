// -*- c++ -*-

#ifndef PYPF_OBJECT_H
#define PYPF_OBJECT_H

#include <string>
#include <map>
#include <mpi.h>
#include "petscfem4py.h"
#include "Options.h"

PYPF_NAMESPACE_BEGIN

class Object {

public:
  virtual ~Object() = 0;
  Object();
  Object(const Object&);


  // reference count
private:
  mutable unsigned int refcnt;
  int get_ref() const { return refcnt;   }
  int inc_ref() const { return ++refcnt; }
  int dec_ref() const { return --refcnt; }
public:
  int getref() const { return get_ref(); }
  int incref() const { return inc_ref(); }
  int decref() const {
    if (get_ref() == 0 || dec_ref() == 0) 
      { delete this; return 0; }
    return get_ref();
  }
public:
  int getRef() const { return get_ref(); }


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
  typedef std::string string;
  bool   hasOption(const string& key) const;
  string getOption(const string& key) const;
  void   setOption(const string& key, const string& value);
  std::map<string,string> getOptions() const;
  void setOptions(const std::map<string,string>&);
  void addOptions(const std::map<string,string>&);

  
  // generic operations
public:
  virtual void setUp()      { }
  virtual void clear()      { }
  virtual void view() const { }
};


PYPF_NAMESPACE_END

#endif // PYPF_OBJECT_H
