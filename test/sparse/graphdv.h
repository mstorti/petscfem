// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: graphdv.h,v 1.2 2002/07/21 16:00:36 mstorti Exp $
#ifndef GRAPHDV_H
#define GRAPHDV_H

#include "graphs.h"
#include "dvector.h"

class int_pair { 
public: 
  int i,j; 
  int_pair(int ii=0,int jj=0) { i=ii; j=jj; }
  bool operator<(int_pair q) const {
    if (i!=q.i) return i<q.i;
    else return j<q.j;
  }
  bool operator==(int_pair q) const {
    return i==q.i && j==q.j;
  }
};

class graphdv : public graph {
private:
  dvector<int_pair> da;
  int ordered;
  int MAX_INIT;
  // will resort if size passes max
  int max;
  // flag whether new elements have been inserted or not
  int modif;
  void resync() { 
    if (modif) {
      // printf("resyncing at size: %d\n",da.size());
      da.sort();
      da.shrink(da.remove_unique());
      modif = 0;
      // printf("end resync.\n",da.size());
    }
  }
  // T* at(int q) { return (T*) da_ref(da,q); }
  void print2dv();
 public:
  graphdv() : da (10000) { 
    ordered = 0; MAX_INIT = 1000; 
    max=MAX_INIT; modif=0; 
  }
  ~graphdv() { da.clear(); }
  void add(int j,int k) {
    int_pair p(j,k);
    if (!da.find(p,0,ordered)) {
      da.push(p); modif=1;      
      int ds = da.size();
      if (ds > max) {
	// launch automatic resync
	resync();
	int new_max = 2 * da.size();
	if (new_max > max) {
	  max = new_max;
	  // printf("new size %d\n",max);
	}
      }
    }  // else printf("already in set...\n");
  }
  void print(const char *s=NULL) { 
    resync(); 
    if (s) printf("%s",s);
    int ss = da.size();
    for (int j=0; j<ss; j++) {
      int_pair &q = da.ref(j);
      printf("(%d,%d) ",q.i,q.j);
    }
  }
  int size() { resync(); return da.size(); }
  void clear() { da.clear(); modif=0; ordered=0; max = MAX_INIT; }
};
#endif
