#include <cstdlib>
#include <cstring>
#include <cstdio>
#include "./hasher.h"

Hasher::Hasher() {
  memset(&xsubi,'\0',
	 3*sizeof(unsigned short int));
  memset(&buffer,'\0',
	 sizeof(struct drand48_data));
}

void Hasher::reset() {
  memset(&xsubi,'\0',
	 3*sizeof(unsigned short int));
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

long int Hasher::val() {
  return result;
}

SumHasher::SumHasher() {
  result = 0 ;
  memset(&buffer,'\0',
	 sizeof(struct drand48_data));
  modulus = 1403711496;
}

void SumHasher::reset() {
  result = 0 ;
}

void SumHasher::hash(int w) {
  long int val;
  memset(&xsubi,'\0',
	 3*sizeof(unsigned short int));
  memcpy(&xsubi,&w,sizeof(int));
  for (int j=0; j<10; j++) 
    nrand48_r(xsubi,&buffer,&val);
  result += val % modulus;
  printf("hash(%d) -> %d, cum %d\n",w,val, result);
}

void SumHasher::hash(int *w,int n) {
  for (int j=0; j<n; j++) {
    hash(w[j]);
  }
}

long int SumHasher::val() {
  return result;
}
