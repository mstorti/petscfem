#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include "./hasher.h"

using namespace std;

void test_hasher(BaseHasher *hash,const char *name) {
  int N=20, M=14, NN=1000000;
  set<int> shash;
  srand(0x234edee4);
  for (int j=0; j<NN; j++) {
    hash->reset();
    for (int k=0; k<N; k++) {
      int w = rand() % M;
      hash->hash(w);
    }
    int s = hash->val();
    shash.insert(s);
  }
  printf("tries %d\n",NN);
  double perf_coll = double(NN)*double(NN)/pow(2.0,32.0)/2.0;
  printf("nbr of collisions: perfect %f, hasher %s: %d\n",
	 perf_coll,name,NN-shash.size());
}

int main() {

#define TEST_HASHER(name)			\
  name name##_h;				\
  test_hasher(&name##_h,#name);

  // TEST_HASHER(Hasher);
  // TEST_HASHER(MD5Hasher);
  // TEST_HASHER(FastHasher);

  TEST_HASHER(SumHasher);
  TEST_HASHER(MD5SumHasher);
  TEST_HASHER(FastSumHasher);
}
