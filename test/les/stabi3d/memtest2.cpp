//__INSERT_LICENSE__
// $Id: memtest2.cpp,v 1.6 2004/02/25 19:27:53 mstorti Exp $
#define _GNU_SOURCE
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include <ctype.h>

using namespace std;

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

  FILE *fid = fopen("/proc/meminfo","r");
  char *line = NULL;
  size_t N=0;
  size_t mem;
  while (1) {
    assert(getline(&line,&N,fid)!=-1);
    if (sscanf(line,"MemFree: %d kB",&mem)) break;
  }
  fclose(fid);
  printf("free: %d\n",mem);
  assert(block_size_m>0.);
  block_size = int(block_size_m*pow(2.,18.));
  nblocks = int(double(mem)/pow(2.,10.)/block_size_m);
  printf("nblocks: %d\n",nblocks);

  MAX = int(double(RAND_MAX)/3.0);
  vector<int *> buffers;
  vector<int> check_sums;
  while (1) {
    int *p = (int *)malloc(sizeof(int)*block_size);
    if(!p) break;
    rand_fill(p,block_size);
    buffers.push_back(p);
    printf("%.2f MB alloc\n",
	   double(buffers.size())
	   *double(block_size)*double(sizeof(int))
	   /pow(2.,20.));
    if (buffers.size() >= nblocks) break;
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
	    failed++;
	  }
	}
      }
    }
    checked++;
    if (checked % 100) printf("checked %d, failed %d\n",checked,failed);
  }
  for (int j=0; j<nblocks; j++) {
    free(buffers[j]);
    buffers[j] = NULL;
  }
}
