// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id merge-with-petsc-233-55-g52bd457 Fri Oct 26 13:57:07 2007 -0300$
#ifndef STAT_H
#define STAT_H

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
  void print_stat(const char * s= NULL) {
    if (s) PetscPrintf(PETSCFEM_COMM_WORLD,
		       "Event %s ------------------------\n",s);
    if (initialized) {
      PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD, 
			      "[%d] total: %g, max: %g, min: "
			      "%g, avrg: %g, count: %d\n",
			      MY_RANK, total(), max(), min(), 
			      avrg(), n());
    } else {
      PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,"[not initialized]\n");
    }
    PetscSynchronizedFlush(PETSCFEM_COMM_WORLD);
  }
};

#endif
