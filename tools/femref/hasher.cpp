#include <cstdlib>
#include <cstring>
#include "./hasher.h"

Hasher::Hasher() {
  memset(&xsubi,'\0',
	 3*sizeof(unsigned short int));
  memset(&buffer,'\0',
	 sizeof(struct drand48_data));
}

void Hasher::reset() {
  memset(&buffer,'\0',
	 sizeof(struct drand48_data));
}

void Hasher::hash(int w) {
  int val;
  memcpy(&val,&xsubi,sizeof(int));
  val += w;
  memcpy(&xsubi,&val,sizeof(int));
  nrand48_r(xsubi,&buffer,&result);
}

void Hasher::hash(int *w,int n) {
  for (int j=0; j<n; j++) hash(w[j]);
}

long int Hasher::hash_val() {
  return result;
}
