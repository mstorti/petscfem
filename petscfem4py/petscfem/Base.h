// -*- c++ -*-
// $Id: Base.h,v 1.1.2.1 2006/03/06 16:56:04 rodrigop Exp $

#ifndef PYPF_SMARTPTR_H
#define PYPF_SMARTPTR_H

#include "namespace.h"
#include "Error.h"

PYPF_NAMESPACE_BEGIN

template<typename T> 
class SmartPtr
{

protected:
  T*  ptr;
  
protected:
  virtual bool valid() {return ptr != NULL; }

public:
  // type aliases
  typedef SmartPtr<T> Ptr;
  typedef T           Base;

public:
  // construction
  explicit SmartPtr(T* t = 0)    : ptr(t)     { }
  SmartPtr(const SmartPtr<T>& p) : ptr(p.ptr) { }
  //destruction
  virtual ~SmartPtr()         { ptr = 0; }
  // assignation
  SmartPtr& operator=(const SmartPtr& sp) { ptr = sp.ptr; return *this; }
  SmartPtr& operator=(const T*& tp)       { ptr = tp;     return *this; }
  // equality
  bool operator==(const SmartPtr& sp) const { return ptr == sp.ptr; }
  bool operator==(const T* tp)        const { return ptr == tp;     }
  bool operator!=(const SmartPtr& sp) const { return ptr != sp.ptr; }
  bool operator!=(const T* tp)        const { return ptr != tp;     }
  // pointer interface
  T* operator->()             { if (!ptr) throw Error("null object"); return  ptr; }
  T& operator*()              { if (!ptr) throw Error("null object"); return *ptr; }
  const T* operator->() const { if (!ptr) throw Error("null object"); return  ptr; }
  const T& operator*()  const { if (!ptr) throw Error("null object"); return *ptr; }
  // automatic conversions
  operator bool()      { return ptr != NULL; }
  operator T&()        { return *ptr;      }
  operator T*&()       { return ptr;       }
  operator const T&()  { return *ptr;      }
  operator const T*&() { return *ptr;      }
  //operator T*()        { return ptr;       }
  //operator const T*()  { return *ptr;      }
  //operator const T*&() { return *ptr;      }
  //operator const T&()  { return *ptr;      }
  //operator T*()   const { return ptr;        }
  //operator T&()   const { return *ptr;       }
  //operator bool() const { return bool(ptr);  }

};

PYPF_NAMESPACE_END


#endif // PYPF_SMARTPTR_H

