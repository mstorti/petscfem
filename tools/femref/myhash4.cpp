#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <vector>
#include <set>
#include <algorithm>
#include "./hasher.h"

using namespace std;

int main() {
  // BasicSumHasher hash;
  // SumHasher hash;
  // MD5SumHasher hash;
  // FastSumHasher hash;
#define SET_HASHER(name) 			\
  name hash;					\
  const char hasher[] = #name

  SET_HASHER(FastHasher);
  time_t start, end;
  int N=20, M=10, NN=100, 
    ntime=10000, NBUFF=N*NN;
  vector<int> buffer(NBUFF);
  for (int j=0; j<NBUFF; j++)
    buffer[j] = rand() % M;
  start = time(NULL);
  for (int l=0; l<ntime; l++) {
    for (int j=0; j<NN; j++) {
      int s=0;
      hash.reset();
      for (int k=0; k<N; k++) {
	hash.hash(&buffer[k*N],N);
      }
      int shv = hash.val();
    }
  }
  double elaps = difftime(time(NULL),start);
  printf("Stats for %s hash class:\n" 
	 "%dx%dx%d evals, elapsed %f, "
	 "rate %g secs/M-int-hash\n",
	 hasher,N,NN,ntime,elaps,
	 elaps/double(NBUFF)/double(ntime)*1.0e6);
}
