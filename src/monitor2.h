// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: monitor2.h,v 1.2 2003/07/08 22:35:52 mstorti Exp $
#ifndef PETSCFEM_MONITOR2_H
#define PETSCFEM_MONITOR2_H

#include <cstdio>
#include <petsc.h>
#include <src/monitor.h>

class DefaultMonitor : public Monitor {
private:
  MPI_Comm comm;
  TextHashTable *options;
  int print_internal_loop_conv;
public:
  void init(MPI_Comm comm,TextHashTable *thash);
  void start();
  void step(int n,double rnorm);
  void stop();
};

#if 0
class ShortMonitor : public Monitor {
};
#endif

#endif
