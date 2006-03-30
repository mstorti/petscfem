// -*- c++ -*-
// $Id: RefCounter.h,v 1.1.2.1 2006/03/30 15:18:14 rodrigop Exp $

#ifndef PYPF_REFCOUNTER_H
#define PYPF_REFCOUNTER_H

#include "namespace.h"

PYPF_NAMESPACE_BEGIN

class RefCounter {

private:
  mutable unsigned int refcnt;
  inline int get_ref() const { return refcnt;   }
  inline int inc_ref() const { return ++refcnt; }
  inline int dec_ref() const { return --refcnt; }

public:
  virtual ~RefCounter() = 0;
  RefCounter() : refcnt(0) { }
  RefCounter(const RefCounter&) : refcnt(0) { }

public:
  inline int getref() const { return get_ref(); }
  inline int incref() const { return inc_ref(); }
  inline int decref() const {
    if (get_ref() == 0 || dec_ref() == 0) { 
      delete this; 
      return 0; 
    }
    return get_ref();
  }

};

PYPF_NAMESPACE_END

#endif // PYPF_REFCOUNTER_H
