// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: buffpack.h,v 1.1 2001/08/01 22:44:21 mstorti Exp $
#ifndef BUFFPACK_H
#define BUFFPACK_H

namespace BufferPack {

  template<class T>
  void pack(const T &t,char *& buff) {
    memcpy(buff,&t,sizeof(T));
    buff += sizeof(T);
  }
  
  template<class T>
  void unpack(T &t,const char *& buff) {
    memcpy(&t,buff,sizeof(T));
    buff += sizeof(T);
  }
  
}

#endif
