//__INSERT_LICENSE__
// $Id: memtest2.cpp,v 1.1 2004/02/25 17:43:34 mstorti Exp $
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
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
  int nblocks = 0;
  int report = 200;
  int regenerate = 20;
  FILE *output = stdout;
  char *log = NULL;
  int len;
  double block_size_m = 1.0;
  int block_size;

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
    case 's':
      sscanf(optarg,"%lf",&block_size_m);
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

  assert(block_size_m>0.);
  block_size = int(block_size*pow(2.,17.));
  assert(ram_size>0.);
  assert(nblocks>0);
  assert(report>0);
  assert(regenerate>0);
  fprintf(output,"checking %.2fMB RAM, %d blocks\n",
	  ram_size,nblocks);

  MAX = int(sqrt(double(RAND_MAX)));
  vector<int *> buffers;
  vector<int> check_sums;
  while (1) {
    int *p = (int *)malloc(sizeof(int)*block_size);
    if(!p) break;
    buffers.push_back(p);
  }
  nblocks = buffers.size();
  while (1) {
    rand_fill(buffers[0],block_size);
    for (int j=1; j<buffers.size(); j++) {
      memcpy(buffers[j],buffers[0],sizeof(int)*block_size);
    }
    // Reordena los buffers;
    random_shuffle(buffers.begin(),buffers.end());
    for (int k=0; k<block_size; k++) {
      int v = buffers[0][k];
      for (int j=1; j<nblocks; j++) {
	if (buffers[j][k]!=v) {
	  if (j<4) break; // may be error in reference
	  else {
	    int w = buffers[j][k];
	    printf("error at %p, old %d, new %d, diff %d\n",
		   &buffers[j][k],v,w,v-w);
	  }
	}
      }
    }
  }
  for (int j=0; j<nblocks; j++) {
    free(buffers[j]);
    buffers[j] = NULL;
  }
}
