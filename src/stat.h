// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: stat.h,v 1.1.120.1 2007/02/24 02:56:42 mstorti Exp $
#ifndef STAT_H
#define STAT_H

extern int MY_RANK,SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class Stat {
public:
  double vmin,vmax,sum;
  int count,initialized;
  Stat() : vmin(NAN), vmax(NAN), sum(NAN), 
           count(0), initialized(0) { }
  void reset() { initialized=0; }
  void add(double val) {
    if (!initialized) {
      initialized = 1;
      vmin = val;
      vmax = val;
      count = 1;
      sum = val;
    } else {
      if (val<vmin) vmin = val;
      if (val>vmax) vmax = val;
      sum += val;
      count++;
    }
  }
  double avrg() { 
    assert(initialized);
    return sum/double(count);
  }
  double min() { 
    assert(initialized);
    return vmin;
  }
  double max() { 
    assert(initialized);
    return vmax;
  }
  int n() {
    assert(initialized);
    return count;
  }
  double total() {
    assert(initialized);
    return sum;
  }
  void print_stat(char * s= NULL) {
    if (s) PetscPrintf(PETSC_COMM_WORLD,
		       "Event %s ------------------------\n",s);
    if (initialized) {
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, 
			      "[%d] total: %g, max: %g, min: "
			      "%g, avrg: %g, count: %d\n",
			      MY_RANK, total(), max(), min(), 
			      avrg(), n());
    } else {
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[not initialized]\n");
    }
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
  }
};

#endif
