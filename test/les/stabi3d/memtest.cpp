//__INSERT_LICENSE__
// $Id: memtest.cpp,v 1.2 2004/02/19 22:11:10 mstorti Exp $
#include <cstdlib>
#include <cstdio>
#include <cmath>

int MAX;

double drand() {
  return double(rand())/double(RAND_MAX);
}

int irand(int M) {
  return int(double(M)*drand());
}

int sum(int *buff,int size) {
  int check = 0;
  for (int j=0; j<size; j++) {
    int k = 2*irand(MAX)-MAX;
    check += k;
    buff[j] = k;
  }
  return check;
}

int main(int argc,char **argv) {
  const double RAM_SIZE = 50;
  const int NBLOCKS = 100;
  MAX = int(sqrt(double(RAND_MAX)));
  int block_size = int(RAM_SIZE*1.0e6/4.0/double(NBLOCKS));
  int size = block_size*NBLOCKS;
  int *buff = new int[size];
  int check = sum(buff,size);
  printf("check: %d\n",check);
  int it=0;
  int bi=0, bj=0;
  while (1) {
    if (bi!=bj) {
      // Intercambia bloques `bi' y `bj'
      int ci = bi*block_size;
      int cj = bj*block_size;
      for (int j=0; j < block_size; j++) {
	int x = buff[ci+j];
	buff[ci+j] = buff[cj+j];
	buff[cj+j] = x;
      }
      int ncheck = sum(buff,size);
      printf("ncheck: %d\n",check);
      if (++it>10) break;
    }
    bi++;
    if (bi>=NBLOCKS) { bi=0; bj++; }
    if (bj>=NBLOCKS) bj=0;
  }
  delete[] buff;
}
