// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: syncbuff.h,v 1.2 2004/01/11 15:01:03 mstorti Exp $
#include <list>
#include <iostream>
#include <src/distcont.h>
#include <src/distcont2.h>
#include <src/debug.h>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <src/autostr.h>

extern int SIZE, MY_RANK;

using namespace std;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<typename T>
class TrivialPartitioner {
public:
  void processor(const T &t,int &nproc,int *plist) {
    nproc=1;
    plist[0] = 0;
  }
};

#define SB(T) DistCont<list<T>,T,TrivialPartitioner<T> >

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<typename T>
class SyncBuffer : public DistCont<list<T>,T,TrivialPartitioner<T> > {
  TrivialPartitioner<T> part;
public:
  SyncBuffer() : 
    DistCont<list<T>,T,TrivialPartitioner<T> >() {}
  void print() { 
    scatter(); 
    sort();
    if (!MY_RANK) {
      typename list<T>::iterator q;
      for (q=begin(); q!=end(); q++) { q->print(); }
    }
  }
  void check_pack();
};


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define SYNC_BUFFER_FUNCTIONS(T)		\
						\
int SB(T)					\
::size_of_pack(const T &t) const {		\
  return t.size_of_pack();			\
}						\
						\
void SB(T)					\
::pack(const T &t, char *&buff) const {		\
  t.pack(buff);					\
}						\
						\
void SB(T)					\
::unpack(T &t,const char *& buff) {		\
  t.unpack(buff); }				\
						\
void SB(T)					\
::combine(const T &t) {				\
push_back(t); };

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<typename T>
void SyncBuffer<T>::check_pack() {
  print();
  int bufsize = 0;
  typename list<T>::iterator q;
  for (q=begin(); q!=end(); q++) 
    bufsize += q->size_of_pack();
  char *buff = new char[bufsize], *p=buff, *pe = buff+bufsize;
  printf("%d elems in buffer\n",list<T>::size());
  for (q =begin(); q!=end(); q++) q->pack(p);
  assert(p == pe);
  clear();
  const char *pp = buff;
  int nelems=0;
  while (pp<pe) {
    nelems++;
    push_back();
    back().unpack(pp);
  }
  assert(pp==pe);
  print();
  printf("%d elems recovered from buffer\n",list<T>::size());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class KeyedLine {
private:
  int key;
  char *line;
public:
  static FILE *output;
  KeyedLine() : key(0), line(NULL) {}
  KeyedLine(const KeyedLine &kl);
  KeyedLine(int key,const AutoString &as);
  ~KeyedLine() { if (line) delete[] line; }
  friend int operator<(const KeyedLine& left, const KeyedLine& right);
  int size_of_pack() const;
  void pack(char *&buff) const;
  void unpack(const char *& buff);
  void print();
};

typedef SyncBuffer<KeyedLine> KeyedOutputBuffer;
