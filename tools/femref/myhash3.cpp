#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include "./hasher.h"

using namespace std;

int count_unique(vector<int> &v) {
  sort(v.begin(),v.end());
  vector<int>::iterator
    last_diff = unique(v.begin(),v.end());
  return last_diff-v.begin();
}

void print_vec(const vector<int> &v) {
  printf("#(");
  for (int k=0; k<v.size(); k++)
    printf("%d ",v[k]);
  printf(")\n");
}

int main() {

  MD5SumHasher hash;
  FastSumHasher phasher;
  int N=20, M=10, NN=10000;
  vector<int>  stat(M);
  set<int> shash, phash;
  for (int j=0; j<NN; j++) {
    hash.reset();
    phasher.reset();
    stat.clear(); stat.resize(M,0);
    for (int k=0; k<N; k++) {
      int w = rand() % M;
      stat[w] += 1;
      hash.hash(w);
      phasher.hash(w);
    }
    int shv = hash.val();
    int s = phasher.val();
    shash.insert(shv);
    phash.insert(s);
  }
  printf("tries %d\n",NN);
  double perf_coll = double(NN)*double(NN)/pow(2.0,32.0);
  printf("nbr of collisions: perfect %f.\n"
	 "For pow-hash-sum %d, md5-hash-sum %d, diff %d\n",
	 perf_coll,NN-phash.size(), NN-shash.size(),
	 phash.size()-shash.size());

}
