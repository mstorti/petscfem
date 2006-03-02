// -*- c++ -*-
// $Id: SmartPtr.h,v 1.1.2.2 2006/03/02 21:37:12 rodrigop Exp $

#ifndef PYPF_SMARTPTR_H
#define PYPF_SMARTPTR_H

#include "namespace.h"
#include "Error.h"

PYPF_NAMESPACE_BEGIN

template<typename T> 
class SmartPtr
{

protected:
  T *ptr;
  
public:
  // type aliases
  typedef SmartPtr<T> Ptr;

public:
  // construction
  explicit SmartPtr(T* t = 0)    : ptr(t)     { }
  SmartPtr(const SmartPtr<T>& p) : ptr(p.ptr) { }
  //destruction
  virtual ~SmartPtr()         { ptr = 0; }
  // assignation
  SmartPtr& operator=(const SmartPtr & sp) { ptr = sp.ptr; return *this; }
  SmartPtr& operator=(T* tp)               { ptr = tp;     return *this; }
  // pointer interface
  T* operator->()             { return  ptr; }
  T& operator*()              { return *ptr; }
  const T* operator->() const { return  ptr; }
  const T& operator*()  const { return *ptr; }
  // automatic conversions
  operator T*()        { return ptr;  }
  operator T*()  const { return ptr;  }
  operator T&()        { return *ptr; }
  operator T&()  const { return *ptr; }

};

PYPF_NAMESPACE_END

#define PYPF_CLASS(NAME) \
class NAME : public SmartPtr< ::NAME >

#define PYPF_CTOR(NAME) \
public: \
NAME(::NAME* p) : Ptr(p) \
{ if (p==NULL) throw Error("null pointer to "#NAME); } \
private:


#endif // PYPF_SMARTPTR_H

