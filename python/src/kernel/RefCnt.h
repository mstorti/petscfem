// $Id$

#ifndef PF4PY_REFCNT_H
#define PF4PY_REFCNT_H

#include "namespace.h"

PF4PY_NAMESPACE_BEGIN

class RefCnt {
private:
  static unsigned long _total;
public:
  static unsigned long allrefs() { return _total;}
private:
  mutable unsigned long _refcnt;
  inline unsigned long get_ref() const { return _refcnt;   }
  inline unsigned long inc_ref() const { ++_total; return ++_refcnt; }
  inline unsigned long dec_ref() const { --_total; return --_refcnt; }
public:
  virtual ~RefCnt() { };
  inline RefCnt() : _refcnt(0) { }
  inline RefCnt(const RefCnt&) : _refcnt(0) { }
  //inline RefCnt& operator=(const RefCnt&) { }
public:
  inline unsigned long getref() const { return get_ref(); }
  inline unsigned long incref() const { return inc_ref(); }
  inline unsigned long decref() const {
    if (get_ref() == 0 || dec_ref() == 0) 
      { delete this; return 0; }
    return get_ref();
  }
};
PF4PY_NAMESPACE_END


PF4PY_NAMESPACE_BEGIN
template<typename T>
class RefVal {
public:
  T* ptr;
  inline ~RefVal() { if (this->ptr) this->ptr->decref(); }
  inline RefVal() : ptr(0)  { }
  inline RefVal(const RefVal& r) : ptr(r.ptr)
  { if (this->ptr) this->ptr->incref(); }
  inline RefVal(const T* t) : ptr(const_cast<T*>(t))
  { if (this->ptr) this->ptr->incref(); }
  inline RefVal(const T& t) : ptr(&const_cast<T&>(t))
  { if (this->ptr) this->ptr->incref(); }
  // assignment
  inline RefVal& operator=(const T* t) {
    if (this->ptr == t) return *this;
    if (this->ptr) this->ptr->decref();
    this->ptr = const_cast<T*>(t);
    if (this->ptr) this->ptr->incref();
    return *this;
  }
  inline RefVal& operator=(const T& t)
  { return this->operator=(&t); }
  inline RefVal& operator=(const RefVal& r)
  { return this->operator=(r->ptr); }
  // equality
  bool operator==(const RefVal& r) const { return this->ptr == r.ptr; }
  bool operator!=(const RefVal& r) const { return this->ptr != r.ptr; }
  bool operator==(const T* t)      const { return this->ptr == t;     }
  bool operator!=(const T* t)      const { return this->ptr != t;     }
  bool operator==(const T& t)      const { return this->ptr == &t;    }
  bool operator!=(const T& t)      const { return this->ptr != &t;    }
  // negation
  inline bool operator!() const { return !this->ptr; }
  // member access
  inline       T* operator->()       { return  this->ptr; }
  inline const T* operator->() const { return  this->ptr; }
  inline       T& operator*()        { return *this->ptr; }
  inline const T& operator*()  const { return *this->ptr; }
  // casting
  inline operator T*() const { return this->ptr; }
  inline operator T&() const { return const_cast<T&>(*this->ptr); }
};
PF4PY_NAMESPACE_END


#include <vector>
PF4PY_NAMESPACE_BEGIN
template<typename T>
class RefVec : private std::vector<T*> {
private:
  RefVec& operator=(const RefVec&);
private:
  typedef typename std::vector<T*>::iterator iter;
  typedef typename std::vector<T*>::const_iterator const_iter;
  inline void incref_all() {
    iter p = this->begin();
    iter e = this->end();
    while (p != e) { const T* t = *p++; t->incref(); }
  }
  inline void decref_all() {
    iter p = this->begin();
    iter e = this->end();
    while (p != e) { const T* t = *p++; t->decref(); }
  }
public:
  using std::vector<T*>::const_iterator;
  using std::vector<T*>::begin;
  using std::vector<T*>::end;
  using std::vector<T*>::size;
  using std::vector<T*>::empty;
public:
  ~RefVec() { this->decref_all(); }
  RefVec() : std::vector<T*>() { }
  RefVec(const RefVec& rv) : std::vector<T*>(rv)
  { this->incref_all(); }
  RefVec(const std::vector<T*>& rv) : std::vector<T*>(rv)
  { this->incref_all(); }
  inline bool operator!() const { return this->empty(); }
  void clear() 
  { this->decref_all(); std::vector<T*>::clear(); }
  T* get(std::size_t k) const {
    if (k < 0 or k >= this->size()) return 0;
    return this->operator[](k);
  }
  bool set(std::size_t k, T* t) {
    if (k < 0 or k >= this->size()) return false;
    if(t == 0) return false;
    t->incref();
    this->operator[](k)->decref();
    this->operator[](k) = t;
    return true;
  }
  bool del(std::size_t k) {
    if (k < 0 or k >= this->size()) return false;
    this->erase(this->begin() + k);
    return true;
  }
  bool has(std::size_t k) const
  { return k > 0 or k < this->size(); }
  bool add(T* t) {
    if (t == 0) return false;
    t->incref();
    this->push_back(t);
    return true;
  }
};
PF4PY_NAMESPACE_END


#include <set>
PF4PY_NAMESPACE_BEGIN
template<typename T>
class RefSet : private std::set<T*> {
private:
  typedef typename std::set<T*>::iterator iter;
  typedef typename std::set<T*>::const_iterator const_iter;
  inline void incref_all() {
    iter s = this->begin();
    iter e = this->end();
    while (s != e) {const T* t = *s++; t->incref();}
  }
  inline void decref_all() {
    iter s = this->begin();
    iter e = this->end();
    while (s != e) {const T* t = *s++; t->decref();}
  }
public:
  using std::set<T*>::const_iterator;
  using std::set<T*>::begin;
  using std::set<T*>::end;
  using std::set<T*>::size;
  using std::set<T*>::empty;
  ~RefSet()  { this->decref_all(); }
  RefSet() : std::set<T*>() { }
  RefSet(const RefSet& rs) : std::set<T*>(rs) 
  { this->incref_all(); }
  RefSet(const std::set<T*>& rs) : std::set<T*>(rs)
  { this->incref_all(); }
  inline bool operator!() const { return this->empty(); }
  void clear() {
    this->decref_all();
    std::set<T*>::clear();
  }
  bool add(T* t) {
    if (t == 0) return false;
    std::pair<iter,bool> p = std::set<T*>::insert(t);
    if (p.second) t->incref();
    return true;
  }
  bool del(T* t) {
    iter s = this->find(t);
    if (s != this->end()) {
      const T* t = *s; t->decref();
      this->erase(s);
      return true;
    }
    return false;
  }
  bool has(T* t) const {
    return this->find(t) != this->end();
  }
};
PF4PY_NAMESPACE_END


#include <map>
PF4PY_NAMESPACE_BEGIN
template<typename K, typename T>
class RefMap : private std::map<K,T*> {
private:
  typedef typename std::map<K,T*>::iterator iter;
  typedef typename std::map<K,T*>::const_iterator const_iter;
public:
  using std::map<K,T*>::const_iterator;
  using std::map<K,T*>::begin;
  using std::map<K,T*>::end;
  using std::map<K,T*>::size;
  using std::map<K,T*>::empty;
public:
  ~RefMap() { this->clear(); }
  RefMap() : std::map<K,T*>() { }
  RefMap(const RefMap& rm) : std::map<K,T*>(rm) {
    const_iter m = this->begin();
    const_iter e = this->end();
    while (m != e) (m++)->second->incref();
  }
  RefMap(const std::map<K,T*>& rm) : std::map<K,T*>(rm) {
    const_iter m = this->begin();
    const_iter e = this->end();
    while (m != e) (m++)->second->incref();
  }
  inline bool operator!() const { return this->empty(); }
  void clear() {
    const_iter m = this->begin();
    const_iter e = this->end();
    while (m != e) (m++)->second->decref();
    std::map<K,T*>::clear();
  }
  T* get(const K& k) const {
    const_iter m = this->find(k);
    if (m != this->end()) return m->second;
    return 0;
  }
  void set(const K& k, T* t) {
    iter m = this->find(k);
    if (m != this->end()) {
      if(t == 0) {
	m->second->decref();
	this->erase(m);
      } else {
	t->incref();
	m->second->decref();
	m->second = t;
      }
    } else {
      if(t == 0) return;
      t->incref();
      this->operator[](k) = t;
    }
  }
  bool del(const K& k) {
    iter m = this->find(k);
    if (m != this->end()) {
      m->second->decref();
      this->erase(m);
      return true;
    }
    return false;
  }
  bool has(const K& k) const {
    return this->find(k) != this->end();
  }
};
PF4PY_NAMESPACE_END


#endif // PF4PY_REFCNT_H

// Local Variables:
// mode: C++
// End:
