// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: buffpack.h,v 1.2 2001/08/06 20:37:12 mstorti Exp $
#ifndef BUFFPACK_H
#define BUFFPACK_H

//#define USE_NAMESPACE

#if USE_NAMESPACE
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
