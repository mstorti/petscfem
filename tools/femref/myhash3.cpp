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

int hash_fun(int w) {
  drand48_data buffer;
  unsigned short int xsubi[3];
  long int val;
  memset(&xsubi,'\0',
	 3*sizeof(unsigned short int));
  memset(&buffer,'\0',
	 sizeof(struct drand48_data));
  memcpy(&xsubi,&w,sizeof(int));
  for (int j=0; j<20; j++) 
    nrand48_r(xsubi,&buffer,&val);
  return val % 10000;
}

int pow_hash(int x) {
  const unsigned int 
    c = 24253464,
    n = 4;
  unsigned long int y, tmp;
  y = x;
  for (int j=0; j<n; j++) {
    tmp = y + c;
    y = tmp*tmp;
  }
  return y;
}

int main() {

#if 1
  MD5SumHasher hash;
  int N=20, M=10, NN=1000000;
  vector<int>  stat(M);
  set<int> shash, sum;
  //#define EXACT_CHECK
#ifdef EXACT_CHECK
  set< vector<int> > val_set;
  map<int, set< vector<int> > > hash_collisions;
#endif
  for (int j=0; j<NN; j++) {
    int s=0;
    hash.reset();
    stat.clear(); stat.resize(M,0);
    for (int k=0; k<N; k++) {
      int w = rand() % M;
      stat[w] += 1;
      hash.hash(w);
      s += pow_hash(w);
    }
    int shv = hash.val();
#if 0
    if (val_set.find(stat) == val_set.end()
	&& (sum.find(s) != sum.end() 
	    || shash.find(shv) != shash.end())) {
      printf("collision! for stat "); print_vec(stat);
      printf(" sum %d, hash-sum %d\n",s,shv);
    }
#endif
#ifdef EXACT_CHECK
    hash_collisions[shv].insert(stat);
    val_set.insert(stat);
#endif
    shash.insert(shv);
    sum.insert(s);
  }
#ifdef EXACT_CHECK
  int total = val_set.size();
  printf("tries %d, different %d\n",NN,total);
  printf("nbr of collisions: for sum %d, hash-sum %d\n",
	 total - sum.size(), total - shash.size());
#else
  printf("tries %d\n",NN);
  double perf_coll = double(NN)*double(NN)/pow(2.0,32.0);
  printf("nbr of collisions: perfect %f.\n"
	 "For sum %d, hash-sum %d, diff\n",
	 perf_coll,NN-sum.size(), NN-shash.size(),
	 sum.size()-shash.size());
#endif

#ifdef EXACT_CHECK
  map<int, set< vector<int> > >::iterator 
    r = hash_collisions.begin();
  while (r != hash_collisions.end()) {
    if (r->second.size()>1) {
      printf("collisions for hash %d:\n",r->first);
      const set< vector<int> > & stat_set = r->second;
      set< vector<int> >::const_iterator w = stat_set.begin();
      while (w != stat_set.end()) print_vec(*w++);
    }
    r++;
  }
#endif
#endif
  
#if 0
  set< vector<int> >::iterator 
    q = val_set.begin();
  printf("differemt combos: \n");
  while (q != val_set.end()) 
    print_vec(*q++);
#endif

#if 0
  int v1[] = {1,1,1,1,1,1,1,2,2,2};
  int v2[] = {0,0,1,1,1,2,2,2,2,2};

  MD5SumHasher hasher;
  // SumHasher hasher;

  hasher.reset();
  hasher.hash(v1,10);
  printf("hash-sum for v1: 0x%08x\n",hasher.val());

  hasher.reset();
  hasher.hash(v2,10);
  printf("hash-sum for v2: 0x%08x\n",hasher.val());
#endif

#if 0
  SumHasher hasher;
  hasher.reset();
  hasher.hash(1);
  hasher.hash(2);
  hasher.hash(3);
  hasher.hash(4);
  hasher.hash(5);
  hasher.hash(6);
  hasher.hash(7);
  hasher.hash(1);
  hasher.hash(2);
  hasher.hash(3);
  hasher.hash(4);
  hasher.hash(5);
  hasher.hash(6);
  hasher.hash(7);
#endif

}
