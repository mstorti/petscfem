// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: fm2temp.h,v 1.1 2003/02/26 22:12:26 mstorti Exp $
#ifndef PETSCFEM_FM2TEMP_H
#define PETSCFEM_FM2TEMP_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This are temporary FastMat2 matrices created on demand. 
    For an object #FastMat2Temp tmp#, #tmp(j)# creates a new 
    #FastMat2# object. Afterwards all this temporary are deleted 
    when the #tmp# object is deleted. 
*/ 
class FastMat2Tmp {
private:
  /// Stores pointers to the matrices
  vector<FastMat2 *> store_v;
  /** For efficiency, we keep a pointer to the store are of the vector. 
      Of course, we must be very careful in keeping this pointer in sync
      with the real store if we resize the vector. This is done in the `sync'
      function. */  
  FastMat2 **store;
  /// The size of the vector
  int size;
  /// Resyncs the pointer #store# and the size to the new size. 
  void sync();
public:
  /** Ctor. */
  FastMat2Tmp();
  /** Dtor. Basically calls clear. */
  ~FastMat2Tmp();
  /** Returns a the j-th matrix. Eventually, it is allocated. 
      @param j (input) index in the vector. */ 
  FastMat2 & operator()(int j);
  /** Cleans up temporary resources. */ 
  void clear();
};

#endif
