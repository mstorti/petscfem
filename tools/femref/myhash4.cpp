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
  MD5SumHasher hash;
  time_t start, end;
  int N=20, M=10, NN=1000, 
    ntime=10, NBUFF=N*NN;
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
  printf("%dx%dx%d evals, elapsed %f, "
	 "rate %g secs/M-int-hash\n",
	 N,NN,ntime,elaps,
	 elaps/double(NBUFF)/double(ntime)*1.0e6);
}
