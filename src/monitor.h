// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: monitor.h,v 1.1 2003/07/07 21:15:26 mstorti Exp $
#ifndef PETSCFEM_MONITOR_H
#define PETSCFEM_MONITOR_H

class Monitor {
 public:
  virtual void start() {}
  virtual void step(int n,double rnorm) {}
  virtual void close() {}
};

#endif
