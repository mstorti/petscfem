#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <set>
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
  SumHasher hash;
  int N=10, M=3, NN=1000;
  vector<int> shash(NN), sum(NN), stat(M);
  set< vector<int> > val_set;
  for (int j=0; j<NN; j++) {
    int s=0;
    hash.reset();
    stat.clear(); stat.resize(M,0);
    for (int k=0; k<N; k++) {
      int w = rand() % M;
      stat[w] += 1;
      hash.hash(w);
      s += w;
    }
    shash[j] = hash.val();
    sum[j] = s;
    val_set.insert(stat);
  }
  int total = val_set.size();
  printf("tries %d, different %d\n",NN,total);
  printf("# of collisions: for sum %d, hash-sum %d\n",
	 total - count_unique(sum), 
	 total - count_unique(shash));
#if 1
  set< vector<int> >::iterator 
    q = val_set.begin();
  printf("differemt combos: \n");
  while (q != val_set.end()) 
    print_vec(*q++);
#endif
}
