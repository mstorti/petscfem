//__INSERT_LICENSE__
// $Id: memtest.cpp,v 1.6 2004/02/24 22:27:26 mstorti Exp $
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <ctype.h>

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

  char *opt_w = NULL;
  char c;
  double ram_size = 0.;
  int nblocks = 100;
  int report = 200;
  int regenerate = 20;
  FILE *output = stdout;
  char *log = NULL;
  int len;

  while ((c = getopt(argc, argv, "l:m:b:p:g:h")) != -1) {
    switch (c) {
    case 'h':
      printf(" usage: $ memtst -m <RAM-in-Mb> "
	     "-b <number-of-blocks> -l <log-file> \n"
	     "-p <report-freq>");
      exit(0);
    case 'm':
      sscanf(optarg,"%lf",&ram_size);
      break;
    case 'b':
      sscanf(optarg,"%d",&nblocks);
      break;
    case 'p':
      sscanf(optarg,"%d",&report);
      break;
    case 'g':
      sscanf(optarg,"%d",&regenerate);
      break;
    case 'l':
      len = strlen(optarg)+1;
      log = new char[len];
      memcpy(log,optarg,len);
      output = fopen(log,"a");
      setvbuf(output,NULL,_IOLBF,0);
      break;
    default:
      if (isprint (optopt))
	fprintf (stderr, "Unknown option `-%c'.\n", optopt);
      else
	fprintf (stderr,
		 "Unknown option character `\\x%x'.\n",
		 optopt);
      abort ();
    }
  }

  assert(ram_size>0.);
  assert(nblocks>0);
  assert(report>0);
  assert(regenerate>0);
  fprintf(output,"checking %.2fMB RAM, %d blocks\n",
	  ram_size,nblocks);

  MAX = int(sqrt(double(RAND_MAX)));
  int block_size = int(ram_size*pow(2.,20.)/4.0/double(nblocks));
  vector<int *> buffers(nblocks);
  vector<int> check_sums(nblocks);
  for (int j=0; j<nblocks; j++) {
    buffers[j] = new int[block_size];
    check_sums[j] = rand_fill(buffers[j],block_size);
  }
  while (1) {
    int block = irand(nblocks);
    int *newblock = new int[block_size];
    memcpy(newblock,buffers[block],sizeof(int)*block_size);
    delete[] buffers[block];
    buffers[block] = newblock;
    int newcheck = sum(buffers[block],block_size) % MAX;
    if (check_sums[block]!=newcheck) {
      fprintf(output,
	      "mem error: old check: %d, new check: %d, diff: %d\n",
	      check_sums[block], newcheck, 
	      newcheck-check_sums[block]);
      failed++;
    }
    check_sums[block] = newcheck;
    checked++;
    if (!irand(regenerate)) {
      check_sums[block] = rand_fill(buffers[block],block_size);
    }
    if (!(checked % report)) {
      double prob = double(failed)/(double(checked)*double(block_size));
      fprintf(output,
	      "checked %d, failed %d, (error prob: %g [#failed/#access])\n",
	      checked, failed, prob);
    }
  }
  for (int j=0; j<nblocks; j++) {
    delete[] buffers[j];
    buffers[j] = NULL;
  }
}
