#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include "./hasher.h"

using namespace std;

int main() {


  MD5Hasher hash;
  FastHasher fhasher;
  Hasher hasher;
  int N=20, M=10, NN=1000000;
  set<int> shash, sfhasher, shasher;
  for (int j=0; j<NN; j++) {
    hash.reset();
    fhasher.reset();
    hasher.reset();
    for (int k=0; k<N; k++) {
      int w = rand() % M;
      hash.hash(w);
      fhasher.hash(w);
      hasher.hash(w);
    }
    int s = hash.val();
    shash.insert(s);

    s = fhasher.val();
    sfhasher.insert(s);

    s = hasher.val();
    shasher.insert(s);
  }
  printf("tries %d\n",NN);
  double perf_coll = double(NN)*double(NN)/pow(2.0,32.0);
  printf("nbr of collisions: perfect %f.\n"
	 "For Md5Hasher %d, FastHasher %d, Hasher %d\n",
	 perf_coll,
	 NN-shash.size(), 
	 NN-sfhasher.size(),
	 NN-shasher.size());
}
