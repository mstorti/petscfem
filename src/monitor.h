// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: monitor.h,v 1.3 2003/07/11 20:47:49 mstorti Exp $
#ifndef PETSCFEM_MONITOR_H
#define PETSCFEM_MONITOR_H

#include <src/texthash.h>

class Monitor {
 public:
  static Monitor *factory(TextHashTable *thash);
  virtual void init(MPI_Comm comm,TextHashTable *thash) {}
  virtual void start() {}
  virtual void step(int n,double rnorm) {}
  virtual void stop() {}
  virtual void close() {}
};

#endif
