// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: debug.h,v 1.2 2001/11/21 19:35:59 mstorti Exp $
#ifndef DEBUG_H
#define DEBUG_H

#include <mpi.h>
#include <string>
#include <map>

/// Puts barriers so that all processes are synchronized
class Debug {
 private:
  /** Flags whether the system should stop or not at each 
      call to `wait' or not 
  */
  map<string,int> active_flags;
  MPI_Comm comm;
  int myrank;
 public:
  int active(const char *s=NULL) const;
  void activate(const char *s=NULL);
  void deactivate(const char *s=NULL);
  void trace(const char *s=NULL);
  Debug(int active_=0,MPI_Comm comm_=MPI_COMM_WORLD);
};

extern Debug debug;

#endif
