// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: timestat.h,v 1.1 2001/11/01 13:56:00 mstorti Exp $
#ifndef TIMESTAT_H
#define TIMESTAT_H

#include <src/fem.h>
//#include <src/readmesh.h>
//#include <src/getprop.h>
#include <src/utils.h>
#include <src/util2.h>
//#include <src/sttfilter.h>
//#include <src/pfmat.h>

#include <src/util2.h>

class TimeStat {
 private:
  double dt,*cumul,*cout;
  int *out,*histo;
  
 public:
  void init();
  void add(double t);
  void print_stat();
  HPChrono chrono;
  double t_min,t_max;
  int nbin;
};

#endif
