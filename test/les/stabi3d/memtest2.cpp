//__INSERT_LICENSE__
// $Id: memtest2.cpp,v 1.7 2004/07/19 11:47:42 mstorti Exp $
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
  double ramp;

  while ((c = getopt(argc, argv, "s:l:m:p:g:h")) != -1) {
    switch (c) {
    case 'h':
      printf(" usage: $ memtest2 -m <%RAM> -l <log-file> -p <report-freq>\n");
      exit(0);
    case 's':
      sscanf(optarg,"%lf",&block_size_m);
      break;
    case 'm':
      sscanf(optarg,"%lf",&ramp);
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
  fprintf(output,"free: %d\n",mem);
  mem = size_t(mem*ramp);
  fprintf(output,"will use: %d\n",mem);
  fclose(fid);
  assert(block_size_m>0.);
  block_size = int(block_size_m*pow(2.,18.));
  fprintf(output,"block_size: %d x sizeof(int)\n",block_size);
  nblocks = int(mem/pow(2.,10.)/block_size_m);
  fprintf(output,"nblocks: %d\n",nblocks);

  MAX = int(double(RAND_MAX)/3.0);
  vector<int *> buffers;
  vector<int> check_sums;
  while (1) {
    int *p = (int *)malloc(sizeof(int)*block_size);
    if(!p) break;
    rand_fill(p,block_size);
    buffers.push_back(p);
    fprintf(output,"%.2f MB alloc\n",
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
	    fprintf(output,"error at %p, old %d, new %d, diff %d\n",
		   &buffers[j][k],v,w,v-w);
	    failed++;
	  }
	}
      }
    }
    checked++;
    if (checked % 100) fprintf(output,"checked %d, failed %d\n",checked,failed);
  }
  for (int j=0; j<nblocks; j++) {
    free(buffers[j]);
    buffers[j] = NULL;
  }
}
