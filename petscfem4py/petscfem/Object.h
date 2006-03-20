// -*- c++ -*-

#ifndef PYPF_OBJECT_H
#define PYPF_OBJECT_H

#include <string>
#include <map>
#include <mpi.h>
#include "petscfem4py.h"

PYPF_NAMESPACE_BEGIN

class Object {

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
  virtual ~Object() = 0;
  Object();
  Object(const Object&);

public:
  int getRef() const { return get_ref(); }
  bool operator==(const Object& obj) const { return this == &obj; }
  bool operator!=(const Object& obj) const { return this != &obj; }


protected:
  MPI_Comm comm;
public:
  MPI_Comm& getComm() const;


protected:
  std::map<std::string,std::string> options;
  virtual OptionTable* get_opt_table() const;
public:
  bool        hasOption(const std::string& key) const;
  std::string getOption(const std::string& key) const;
  void        setOption(const std::string& key,
			const std::string& value);
  std::map<std::string,std::string> getOptions() const;
  void setOptions(const std::map<std::string,std::string>&);

  
public:
  virtual void setUp()      { }
  virtual void clear()      { }
  virtual void view() const { }

};


PYPF_NAMESPACE_END

#endif // PYPF_OBJECT_H
