#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include "./hasher.h"

using namespace std;

int main() {
  Hasher hash;
  for (int N=16; N<24; N++) {
    // printf("N %d\n",N);
    vector<int> vec(N);
    vector<int> hashset;
    for (int j=0; j<N; j++) vec[j]=0;
    int count=0;
    int print=0;
    while (1) {
      hash.reset();
      hash.hash(&vec[0],vec.size());
      int h = hash.hash_val();
      if (print) {
	printf("vec: ");
	for (int j=0; j<N; j++) 
	  printf("%d ",vec[j]);
	printf(", hash: %x\n",h);
      }
      hashset.push_back(h);
      count++;
      int j;
      for (j=0; j<N; j++) {
	if (vec[j]) vec[j]=0;
	else {
	  vec[j]=1;
	  break;
	}
      }
      if (j>=N) break;
    }
    sort(hashset.begin(),hashset.end());
    vector<int>::iterator
      last_diff = unique(hashset.begin(),hashset.end());
    int ndiff = last_diff-hashset.begin();
    printf("N %d, %d bit vectors checked,"
	   " %d different hash values, "
	   "%d collisions\n",
	   N,count,ndiff,count-ndiff);
  }
}
