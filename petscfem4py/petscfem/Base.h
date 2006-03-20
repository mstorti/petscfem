// -*- c++ -*-
// $Id: Base.h,v 1.1.2.2 2006/03/20 16:06:00 rodrigop Exp $

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
  typedef SmartPtr<T> Ptr;
  typedef T           Base;

public:
  //destruction
  virtual ~SmartPtr() { ptr = 0; }
  // construction
  explicit SmartPtr(T* t, bool flag=true) : ptr(t)
  { if (flag && ptr == NULL) throw Error("null pointer"); }
  SmartPtr(const SmartPtr<T>& p) : ptr(p.ptr)
  { if (ptr == NULL) throw Error("null pointer"); }
  // assignation
  SmartPtr& operator=(const T*& tp)  { ptr = tp; return *this; }
  // equality
  //bool operator==(const SmartPtr& sp) const { return ptr == sp.ptr; }
  //bool operator==(const T* tp)        const { return ptr == tp;     }
  //bool operator!=(const SmartPtr& sp) const { return ptr != sp.ptr; }
  //bool operator!=(const T* tp)        const { return ptr != tp;     }
  // pointer interface
  T* operator->()
  { if (!ptr) throw Error("null pointer"); return  ptr;  }
  const T* operator->() const
  { if (!ptr) throw Error("null pointer"); return  ptr;  }
  // automatic conversions
  operator T&()         { return *ptr; }
  //operator T*&()        { return ptr;  }
  //operator T&()  const { return *ptr; }
  operator T* ()             { return ptr;  }
  //operator T*& ()            { return ptr;  }
  operator T* const () const { return ptr;  }

};

PYPF_NAMESPACE_END


#endif // PYPF_SMARTPTR_H

