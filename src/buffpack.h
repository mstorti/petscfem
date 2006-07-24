// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: buffpack.h,v 1.5 2006/07/24 04:28:15 mstorti Exp $
#ifndef BUFFPACK_H
#define BUFFPACK_H

#define BUFFER_PACK BufferPack::pack
#define BUFFER_UNPACK BufferPack::unpack

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
