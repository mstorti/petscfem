//__INSERT_LICENSE__
// $Id: memtest.cpp,v 1.4 2004/02/24 18:16:57 mstorti Exp $
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>

int MAX;

double drand() {
  return double(rand())/double(RAND_MAX);
}

int irand(int M) {
  return int(double(M)*drand());
}

int sum(int *buff,int size) {
  int check = 0;
  for (int j=0; j<size; j++) 
    check = (check + buff[j]) % MAX;
  return check;
}

int rand_fill(int *buff,int block_size) {
  for (int k=0; k<block_size; k++) 
    buff[k] = irand(MAX);
  return sum(buff,block_size) % MAX;
}

int main(int argc,char **argv) {
  int checked=0, failed=0;
  int report = 200;
  const double RAM_SIZE = 480;
  const int NBLOCKS = 100;
  MAX = int(sqrt(double(RAND_MAX)));
  int block_size = int(RAM_SIZE*pow(2.,20.)/4.0/double(NBLOCKS));
  vector<int *> buffers(NBLOCKS);
  vector<int> check_sums(NBLOCKS);
  for (int j=0; j<NBLOCKS; j++) {
    buffers[j] = new int[block_size];
    check_sums[j] = rand_fill(buffers[j],block_size);
  }
  while (1) {
    int block = irand(NBLOCKS);
    int *newblock = new int[block_size];
    memcpy(newblock,buffers[block],sizeof(int)*block_size);
    delete[] buffers[block];
    buffers[block] = newblock;
    int newcheck = sum(buffers[block],block_size) % MAX;
    if (check_sums[block]!=newcheck) {
      printf("mem error: old check: %d, new check: %d, diff: %d\n",
	     check_sums[block], newcheck, newcheck-check_sums[block]);
      failed++;
    }
    check_sums[block] = newcheck;
    checked++;
    if (!irand(20)) {
      // printf("regenera bloque %d\n",block);
      check_sums[block] = rand_fill(buffers[block],block_size);
    }
    if (!(checked % report)) 
      printf("checked %d, failed %d\n",checked, failed);
  }
  for (int j=0; j<NBLOCKS; j++) {
    delete[] buffers[j];
    buffers[j] = NULL;
  }
}
