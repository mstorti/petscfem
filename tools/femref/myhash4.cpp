#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <vector>
#include <set>
#include <algorithm>
#include "./hasher.h"

using namespace std;

void test_hasher(BaseHasher *hash,
		 const char *name) {
  time_t start, end;
  int N=20, M=1000, NN=100, 
    ntime=10000, NBUFF=N*NN;
  vector<int> buffer(NBUFF);
  for (int j=0; j<NBUFF; j++)
    buffer[j] = rand() % M;
  start = time(NULL);
  for (int l=0; l<ntime; l++) {
    for (int j=0; j<NN; j++) {
      int s=0;
      hash->reset();
      for (int k=0; k<N; k++) {
	hash->hash(&buffer[k*N],N);
      }
      int shv = hash->val();
    }
  }
  double elaps = difftime(time(NULL),start);
  printf("Stats for %s hash class:\n" 
	 "%dx%dx%d evals, elapsed %f, "
	 "rate %g secs/M-int-hash\n",
	 name,N,NN,ntime,elaps,
	 elaps/double(NBUFF)/double(ntime)*1.0e6);
}

int main() {
#define DEF_HASHER(type) \
  type type##_h

  DEF_HASHER(Hasher);
  DEF_HASHER(SumHasher);
  DEF_HASHER(MD5Hasher);
  DEF_HASHER(MD5SumHasher);
  DEF_HASHER(FastHasher);
  DEF_HASHER(FastSumHasher);
  DEF_HASHER(BJHasher);

  const char *hasher;

#define TEST_HASHER(name)				\
  test_hasher(&name##_h,#name)

  TEST_HASHER(Hasher);
  TEST_HASHER(MD5Hasher);
  TEST_HASHER(FastHasher);

  TEST_HASHER(SumHasher);
  // TEST_HASHER(MD5SumHasher);
  TEST_HASHER(FastSumHasher);
  TEST_HASHER(BJHasher);

}
