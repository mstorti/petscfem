#include <cstdlib>
#include <cstring>
#include <vector>
#include <set>

using namespace std;

int myhash(int *w,int n) {
  drand48_data buffer;
  unsigned short int xsubi[3];
  memset(&buffer,'\0',sizeof(struct drand48_data));
  memset(&xsubi,'\0',3*sizeof(unsigned short int));
  long int result;
  int val;
  for (int j=0; j<n; j++) {
    memcpy(&val,&xsubi,sizeof(int));
    val += w[j];
    memcpy(&xsubi,&val,sizeof(int));
    nrand48_r(xsubi,&buffer,&result);
  }
  return result;
}

int main() {
  for (int N=16; N<24; N++) {
    // printf("N %d\n",N);
    vector<int> vec(N);
    set<int> hashset;
    for (int j=0; j<N; j++) vec[j]=0;
    int count=0;
    int print=0;
    while (1) {
      int h = myhash(&vec[0],N);
      if (print) {
	printf("vec: ");
	for (int j=0; j<N; j++) 
	  printf("%d ",vec[j]);
	printf(", hash: %x\n",h);
      }
      hashset.insert(h);
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
    printf("N %d, %d bit vectors checked,"
	   " %d different hash values, "
	   "%d collisions\n",
	   N,count,hashset.size(),count-hashset.size());
  }
}
