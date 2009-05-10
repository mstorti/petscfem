// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id config-1.0.7-74-gc7a455f Mon Oct 29 23:28:05 2007 -0300$
#ifndef PETSCFEM_FM2STATS_H
#define PETSCFEM_FM2STATS_H

class FastMat2Stats {
 public:
  int 
    use_dgemm,
    was_sl,
    was_sl_count,
    was_not_sl_count;
 FastMat2Stats() 
   : use_dgemm(1), 
    was_sl(0),
    was_sl_count(0),
    was_not_sl_count(0) {}
  void report() {
    int total = 
      was_not_sl_count + was_sl_count;
    double ratio = 100.0*double(was_sl_count)/total;
    printf("FM2STATS: sl %d(%.3f%%), not sl %d, total %d\n",
           was_sl_count,ratio,was_not_sl_count,total);
  }
  ~FastMat2Stats() { report(); }
};

extern FastMat2Stats glob_fm2stats;

#endif
