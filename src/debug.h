// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: debug.h,v 1.1 2001/11/20 02:28:59 mstorti Exp $
#ifndef DEBUG_H
#define DEBUG_H

/// Puts barriers so that all processes are synchronized
class Debug {
 private:
  /** Flags whether the system should stop or not at each 
      call to `wait' or not 
  */
  int active_;
  MPI_Comm comm;
 public:
  int active() const { return active_;} 
  void activate() { active_=1;}
  void deactivate() { active_=0;}
  void wait(const char *s=NULL);
  Debug(int active__=0,MPI_Comm comm_=MPI_COMM_WORLD) : 
    active_(active__), comm(comm_) {}
};

extern Debug debug;

#endif
