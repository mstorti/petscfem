#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include "./hasher.h"

using namespace std;

int main() {
  Hasher hash;
  for (int N=3; N<4; N++) {
    vector<int> vec(N);
    int NN = int(pow(2.0,N));
    vector<int> hashset(NN);
    for (int step=0; step<2; step++) {
	// printf("N %d\n",N);
      for (int j=0; j<N; j++) vec[j]=0;
      int count=0;
      int print=0;
      while (1) {
#if 1
	hash.reset();
	hash.hash(&vec[0],vec.size());
	int h = hash.hash_val();
	if (print) {
	  printf("vec: ");
	  for (int j=0; j<N; j++) 
	    printf("%d ",vec[j]);
	  printf(", hash: %x\n",h);
	}
	if (step==0)
	  hashset[count] = h;
	else printf("hash before %d, now %d\n",
		    hashset[count],hash.hash_val());
#else
	hashset[count] = rand();
#endif
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
    }
  }
}
