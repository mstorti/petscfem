// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: buffpack.h,v 1.3 2001/08/07 16:57:57 mstorti Exp $
#ifndef BUFFPACK_H
#define BUFFPACK_H

//#define USE_NAMESPACE

#ifdef USE_NAMESPACE
#define BUFFER_PACK BufferPack::pack
#define BUFFER_UNPACK BufferPack::unpack
#else
#define BUFFER_PACK BufferPack_pack
#define BUFFER_UNPACK BufferPack_unpack
#endif

#ifdef USE_NAMESPACE
namespace BufferPack {
#endif

  template<class T>
  void BUFFER_PACK(const T &t,char *& buff) {
    memcpy(buff,&t,sizeof(T));
    buff += sizeof(T);
  }
  
  template<class T>
  void BUFFER_UNPACK(T &t,const char *& buff) {
    memcpy(&t,buff,sizeof(T));
    buff += sizeof(T);
  }
  
#ifdef USE_NAMESPACE
}
#endif

#endif
