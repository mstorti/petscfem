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
  map<unsigned int,unsigned int> table;
  int npoints, nint;
  double tol;
  int ndim;
  MD5Hasher hasher;
  vector<double> bounds;
  double & bound(int ndim, int which) {
    return bounds[3*ndim+which];
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
  int add(const vector<double> &x) {
    assert(x.size()%ndim==0);
    int np = x.size()/ndim;
    int ok = 0;
    for (int j=0; j<np; j++) {
      hasher.reset();
      for (int k=0; k<ndim; k++) {
	double xx = x[j*ndim+k];
	coords.push_back(xx);
	unsigned int w = 
	  (unsigned int)(nint*(xx-bound(k,0))/bound(k,2));
	// printf(" %x",w);
	hasher.hash(w);
      }
      // printf("\n");
      int h = hasher.val();
//       printf("adding point j %d, (%f,%f,%f), hash %x\n",
// 	     j,x[j*ndim],x[j*ndim+1],x[j*ndim+2],h);
      if (table.find(h)==table.end()) {
	table[h] = npoints++;
	ok ++;
      } 
    }
    return ok;
  }

  int get(const vector<double> &x,
	  vector<double> &ngbrs) {
    int np = x.size()/ndim;
    for (int j=0; j<np; j++) {
      hasher.reset();
      for (int k=0; k<ndim; k++) {
	double xx = x[j*ndim+k];
	int w = int(nint*(xx-bound(j,0))/bound(j,2));
	hasher.hash(w);
      }
      int h = hasher.val();
      if (table.find(h)!=table.end()) {
	int indx = table[h];
	for (int j=0; j<ndim; j++) {
	  ngbrs.push_back(coords[indx*ndim+j]);
	}
      }
    }
  }
};

double drand() {
  return double(rand())/double(RAND_MAX);
}

// Tries to solve the ANN problem (or related)
// through hashing
int main() {
  int N=10000000, ndim=3;			// Number of points to be added
  hashed_coords_t hashed_coords(ndim,1e-10);
  vector<double> coords;
  for (int j=0; j<N*ndim; j++)
    coords.push_back(drand());
  int ok = hashed_coords.add(coords);
  printf("tried %d, ok %d\n",N,ok);
  int bad=0;
  for (int j=0; j<N; j++) {
    vector<double> xtry;
    int q = rand() % 2;
    if (q) {
      // printf("random generation\n");
      for (int j=0; j<ndim; j++) 
	xtry.push_back(drand());
    } else {
      int k = rand()%N;
      // printf("probing point %d\n",k);
      for (int j=0; j<ndim; j++) 
	xtry.push_back(coords[k*ndim+j]);
    }
    vector<double> found;
    hashed_coords.get(xtry,found);
    bad += (found.size()==0 != q);
  }
  printf("total %d, OK %d, bad %d\n",N,N-bad,bad);
}
