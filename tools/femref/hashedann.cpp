#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cassert>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>

#include "./hasher.h"

using namespace std;

class hashed_coords_t {
private:
  vector<double> coords;
  typedef 
  multimap<unsigned int,unsigned int> map_t;
  typedef 
  pair<map_t::const_iterator,map_t::const_iterator> pair_iter_t;
  typedef pair<unsigned int,unsigned int> pair_t;
  map_t table;
  int npoints, nint;
  double tol;
  int ndim;
  MD5Hasher hasher;
  vector<double> bounds;
  double & bound(int ndim, int which) {
    return bounds[3*ndim+which];
  }
  void get_candidates(const vector<int> &intervals,
		      const vector<int> &incs, 
		      vector<int> &candidates) {
    hasher.reset();
    for (int j=0; j<ndim; j++) 
      hasher.hash(intervals[j]+incs[j]);
    int h = hasher.val();
    pair_iter_t p = table.equal_range(h);
    for (map_t::const_iterator q = p.first; q!=p.second; q++)
      candidates.push_back(q->second);
  }
  void get_candidates(const double *x,
		      vector<int> &candidates) {
    vector<int> intervals(ndim);
    for (int k=0; k<ndim; k++) {
      intervals[k] = 
	(unsigned int)(nint*(x[k]-bound(k,0))/bound(k,2));
    }
    vector<int> incs(ndim);
    for (int k=0; k<ndim; k++) incs[k] = -1;
    while (1) {
      int ok=0;
      for (int j=0; j<ndim; j++) {
	if (incs[j]<1) { ok=1; incs[j]++; break; } 
	else { incs[j]=-1; }
      }
      if (!ok) break;
      get_candidates(intervals,incs,candidates);
    }
  }
  void get_candidates(const vector<double> &x,
		      vector<int> &candidates) {
    assert(x.size() % ndim == 0);
    int np = x.size()/ndim;
    candidates.clear();
    for (int j=0; j<np; j++) 
      get_candidates(&x[j*ndim],candidates);
  }
  double dist(int k,const vector<double> &x) {
    double d = 0.0;
    for (int j=0; j<ndim; j++) {
      double dx = coords[k*ndim+j]-x[j];
      d += dx*dx;
    }
    d = sqrt(d);
    return d;
  }
public:
  hashed_coords_t(int ndim_a, double tol_a=NAN)
    : ndim(ndim_a), tol(tol_a), npoints(0) {
    if (isnan(tol)) tol=0.1/UINT_MAX;
    bounds.resize(3*ndim);
    for (int j=0; j<ndim; j++) {
      bound(j,0) = 0.0;
      bound(j,1) = 1.0;
      bound(j,2) = 1.0;
    }
  }
  int add(const vector<double> &x,double tol_a=NAN) {
    if (isnan(tol_a)) tol_a=tol;
    assert(tol_a<=tol);
    assert(x.size() % ndim == 0);
    int np = x.size()/ndim;
    int ok = 0;
    for (int j=0; j<np; j++) {
      vector<int> candidates;
      vector<double> xx(ndim);
      for (int k=0; k<ndim; k++) 
	xx[k] = x[j*ndim+k];
      get_candidates(xx,candidates);
      int ncand = candidates.size();
      for (int j=0; j<ncand; j++) {
	double d = dist(candidates[j],x);
	if (d<tol_a) continue;
      }
      ok++;
      hasher.reset();
      for (int k=0; k<ndim; k++) {
	coords.push_back(xx[k]);
	int interv = 
	  (unsigned int)(nint*(xx[k]-bound(k,0))/bound(k,2));
	hasher.hash(interv);
      }
      int w = hasher.val();
      table.insert(pair_t(w,npoints++));
    }
    return ok;
  }
  int size() { return npoints; }
  int get(const vector<double> &x,
	  vector<double> &ngbrs,double tol_a=NAN) {
    if (isnan(tol_a)) tol_a=tol;
    assert(tol_a<=tol);
    assert(x.size()==ndim);
    vector<int> candidates;
    get_candidates(x,candidates);
    int ncand = candidates.size();
    int OK=0;
    for (int k=0; k<ncand; k++) {
      int c = candidates[k];
      double d = dist(c,x);
      if (d<=tol_a) {
	OK++;
	for (int j=0; j<ndim; j++) 
	  ngbrs.push_back(coords[c*ndim+j]);
      }
    }
    return OK;
  }
  void print() {
    map_t::iterator q = table.begin();
    while (q!=table.end()) {
      printf("%d -> %d\n",q->first,q->second);
      q++;
    }
  }
};

double drand() {
  return double(rand())/double(RAND_MAX);
}

// Tries to solve the ANN problem (or related)
// through hashing
int main() {
  int N=10, ndim=3;			// Number of points to be added
  hashed_coords_t hashed_coords(ndim,1e-10);
  vector<double> coords;
  for (int j=0; j<N*ndim; j++)
    coords.push_back(drand());
  int npoints = hashed_coords.add(coords);
  printf("tried %d, OK %d\n",N,npoints);
  int bad=0;
  hashed_coords.print();
  for (int j=0; j<N; j++) {
    vector<double> xtry;
    int q = rand() % 2;
    if (q) {
      printf("random generation\n");
      for (int j=0; j<ndim; j++) 
	xtry.push_back(drand());
    } else {
      int k = rand()%N;
      printf("probing point %d\n",k);
      for (int j=0; j<ndim; j++) 
	xtry.push_back(coords[k*ndim+j]);
    }
    vector<double> found;
    hashed_coords.get(xtry,found);
    bad += (found.size()==0 != q);
  }
  printf("total %d, OK %d, bad %d\n",N,N-bad,bad);
}
