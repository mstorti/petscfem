//__INSERT_LICENSE__
// $Id: socket2.cpp,v 1.6 2003/02/03 18:35:25 mstorti Exp $
#define _GNU_SOURCE
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cmath>
#include <HDR/sockets.h>

enum comm_mode { SEND, RECV };

int count, inside = 0;
const int CHUNK=1000;

void chomp(char *s) { s[strlen(s)-1] = '\0'; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define SGETLINE_FACTOR 2
// #define SGETLINE_INIT_SIZE 128
#define SGETLINE_INIT_SIZE 8
#define SGETLINE_MAX_SIZE INT_MAX
ssize_t Sgetline(char **lineptr, size_t *N_a,Socket *sock) {
  unsigned int &N = *N_a;	// better readbility
  char * new_line_ptr = NULL, *q, *qe;
  int old_N;
  while (1) {
    old_N = N;
    N = (N ? 2*N : SGETLINE_INIT_SIZE);
#define DEBUG
#ifdef DEBUG
    printf("Allocating %d bytes\n",N);
#endif
    if (N > SGETLINE_MAX_SIZE) return 0;
    new_line_ptr = (char *) malloc(N);
    assert(new_line_ptr);
    if (*lineptr) {
      memcpy(new_line_ptr,*lineptr,old_N);
      free(*lineptr);
    }
    *lineptr = new_line_ptr;
    new_line_ptr = NULL;
    Sgets(*lineptr+old_N,N-old_N,sock);
    qe = *lineptr+N;
    for (q = *lineptr+old_N; q<qe; q++) {
      if (*q == '\0') break;
    }
    if (q<qe) break;
  }
  return strlen(*lineptr)+1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
inline double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int talk(Socket *sock,comm_mode &mode) { 
  static char *line = NULL;
  static size_t N=0,NN=0;
  char *buf=NULL;
  int stop;

  if (mode==SEND) {
    printf("> ");
    getline (&line,&N,stdin);
    Sprintf(sock,line);
    if (!strcmp(line,"OVER\n")) mode = RECV;
    stop = !strcmp(line,"STOP\n");
  } else if (mode==RECV) {
    int computed=0;
    while (Stest(sock)==0) {
      for (int k=0; k<CHUNK; k++) {
	double x,y;
	x = drand();
	y = drand();
	inside += (x*x+y*y<1.);
      }
      computed += CHUNK;
      count += CHUNK;
    }
    double PI = 4.*double(inside)/double(count);
    printf("[Computed %d points, current PI %f, error %g]\n",computed,
	   PI,fabs(PI-M_PI));
    
    Sgetline(&buf,&NN,sock);
    if (!strcmp(buf,"OVER\n")) mode = SEND;
    printf("-- %s",buf);
    stop = !strcmp(buf,"STOP\n");
  } else assert(0);
  return stop;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int main(int argc,char **args) {
  comm_mode mode;
  Socket *SRVR=NULL, *srvr = NULL, *clnt = NULL;
  if(argc>1 && !strcmp(args[1],"-server")) {
    SRVR = Sopen("","s5555");
    assert(SRVR);
    printf("Waiting connection...\n");
    srvr = Saccept(SRVR);
    assert(srvr);
    mode = SEND;
    while (!talk(srvr,mode));
    Sclose(srvr);
    Sclose(SRVR);
  } else {
    clnt = Sopen("","c5555");
    assert(clnt);
    mode = RECV;
    while (!talk(clnt,mode));
    Sclose(clnt);
  }
  exit(EXIT_SUCCESS);
}
