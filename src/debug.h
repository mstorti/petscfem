// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: debug.h,v 1.9 2003/07/02 02:32:47 mstorti Exp $
#ifndef PF_DEBUG_H
#define PF_DEBUG_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <unistd.h>
#include <sys/types.h>

#include <signal.h>

#include <mpi.h>
#include <string>
#include <map>

#include <src/utils.h>
#include <src/util2.h>

/// Puts barriers so that all processes are synchronized
class Debug {
 private:
  /** Flags whether the system should stop or not at each 
      call to `wait' or not 
  */
  map<string,int> active_flags;
  MPI_Comm comm;
  int myrank,size;
  HPChrono chrono;
  static int stop_f;
  static sighandler_t orig_handler;
  char *line;
  size_t N;
  int was_initialized;
  vector<int> flags;
  void release_proc(int proc=0);
  int dummy;
 public:
  static void init();
  static void set_signal(int signal);
  int active(const char *s=NULL) const;
  void activate(const char *s=NULL);
  void deactivate(const char *s=NULL);
  void trace(const char *s=NULL);
  Debug(int active_=0,MPI_Comm comm_=PETSC_COMM_WORLD);
};

extern Debug *GLOBAL_DEBUG;
#endif

