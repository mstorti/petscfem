// -*- c++ -*-
// $Id: SmartPtr.h,v 1.1.2.4 2006/03/30 15:18:14 rodrigop Exp $

#ifndef PYPF_SMARTPTR_H
#define PYPF_SMARTPTR_H

#include "namespace.h"
#include "Error.h"

PYPF_NAMESPACE_BEGIN

template<typename T> 
class SmartPtr
{

private:
  SmartPtr(): ptr(0) { }
  
protected:
  T* ptr;
  
public:
  // type aliases
  typedef SmartPtr<T> Handle;
  typedef T           Base;

public:
  // destruction
  virtual ~SmartPtr() { this->ptr = 0; }
  // construction
  explicit SmartPtr(T* t, bool flag=true) : ptr(t) { 
    if (flag && this->ptr == NULL) 
      throw Error("constructing a handle from a null pointer"); 
  }
  SmartPtr(const SmartPtr<T>& p) : ptr(p.ptr){ 
    if (this->ptr == NULL)
      throw Error("constructing a handle from a null pointer");
  }
  // assignment
  inline SmartPtr& operator=(T* tp) { 
    this->ptr = tp;
    return *this;
  }
  inline SmartPtr& operator=(const T*& tp) { 
    this->ptr = tp;
    return *this;
  }
  inline SmartPtr& operator=(const SmartPtr<T>& sp) {
    this->ptr = sp.ptr;
    return *this;
  }
  // equality/not equal
  bool operator==(const T*& tp)       const { return this->ptr == tp;     }
  bool operator!=(const T*& tp)       const { return this->ptr != tp;     }
  bool operator==(const SmartPtr& sp) const { return this->ptr == sp.ptr; }
  bool operator!=(const SmartPtr& sp) const { return this->ptr != sp.ptr; }
  // negation
  inline bool operator!() {
    return !this->ptr;
  }
  // member access
  inline T* operator->() {
    if (!this->ptr) throw Error("trying to access a null pointer"); 
    return  this->ptr;
  }
  inline const T* operator->() const { 
    if (!this->ptr) throw Error("trying to access a null pointer"); 
    return  ptr;
  }
  // casting
  inline operator       T&()        { return *this->ptr; }
  inline operator const T&() const  { return *this->ptr; }
  inline operator T*&()             { return this->ptr;  }
  inline operator T* const&() const { return this->ptr;  }

};

PYPF_NAMESPACE_END

#endif // PYPF_SMARTPTR_H
