// -*- mode: C++ -*- 
#ifndef WSPACE_H
#define WSPACE_H

#include <set>

template<class T>
class WorkSpace : public set<T *> {
public:
  T *newobj() {T *t = new T; insert(t); return t;}
  void free(T *t) {erase(t);}
  ~WorkSpace();
};

template<class T>
WorkSpace<T>::~WorkSpace() {
  for (set<T *>::iterator j=begin(); j!=end(); j++) {
    delete *j;
  }
}

#endif
