// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: sockbuff.h,v 1.1 2003/06/08 13:10:43 mstorti Exp $
#ifndef PETSCFEM_SOCKBUFF_H
#define PETSCFEM_SOCKBUFF_H
#include <src/dvector.h>

#define BUFF_SIZE 50000
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
class SocketBuffer {
private:
  dvector<T> buff;
  int last;
  Socket *sock;

public:
  SocketBuffer(Socket *sock_a, int size=BUFF_SIZE) 
    : last(0), buff(size), sock(sock_a) {
    buff.resize(size);
  }

  ~SocketBuffer() { buff.clear(); }

  void flush() {
    if (last) Swrite(sock,buff.buff(),last*sizeof(T));
    last = 0;
  }
  
  void put(T t) {
    buff.e(last++) = t;
    if (last>=buff.size()) flush();
  }
};

#endif
