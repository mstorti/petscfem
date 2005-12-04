// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: iisdmatstat.h,v 1.1 2005/12/04 10:03:19 mstorti Exp $
#ifndef PETSCFEM_IISDMATSTAT_H
#define PETSCFEM_IISDMATSTAT_H

struct iisdmat_stat_t {
  // averg, min, max, 
  double local[3], interf[3];
  int count;
  iisdmat_stat_t();
  void reset();
  void report();
  void stat(double local, double interf);
};

extern iisdmat_stat_t iisdmat_stat;

#endif
