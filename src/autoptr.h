// -*- mode: C++ -*-
//__INSERT_LICENSE__
// $Id: autoptr.h,v 1.1 2003/02/16 22:03:18 mstorti Exp $
#ifndef PETSCFEM_AUTOPTR_H
#define PETSCFEM_AUTOPTR_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
class AutoPtr {
private:
  T *array;
  int rank;
  vector<int> shape;
public:
  AutoString();
  ~AutoString();
  T *store();
  T& operator[]();
  T& operator[](int j,...);
  const T& operator[]() const;
  T& operator[](int j,...) const;
  int size();
  void reshape(int rank,...);
  void resize(int rank,...);
  void clear();
};

#endif
