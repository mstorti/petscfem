//__INSERT_LICENSE__
// $Id: socket2.cpp,v 1.9 2003/02/08 23:12:57 mstorti Exp $
#define _GNU_SOURCE
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cmath>
#include <HDR/sockets.h>
#include <pthread.h>

enum comm_mode { SEND, RECV };

int count, inside = 0;
const int CHUNK=1000;

void chomp(char *s) { s[strlen(s)-1] = '\0'; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void pprint(char *buf,int N) {
  char *q = buf;
  char *qe = q+N;
  for (q = buf; q<qe; q++) {
    if (*q == '\n') {
      putc('\\',stdout);
      putc('n',stdout);
    } else if (*q == '\0') {
      putc('\\',stdout);
      putc('0',stdout);
    } else putc(*q,stdout);
  }
  putc('\n',stdout);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define SGETLINE_FACTOR 2
// #define SGETLINE_INIT_SIZE 128
#define SGETLINE_INIT_SIZE 16
#define SGETLINE_MAX_SIZE INT_MAX
ssize_t Sgetline(char **lineptr, size_t *N_a,Socket *sock) {
  unsigned int &N = *N_a;	// better readbility
  char * new_line_ptr = NULL, *q, *q0, *qe;
  int read_so_far=0;
  while (1) {
    if (N>0) {
      q0 = *lineptr+read_so_far;
      qe = *lineptr+N;
      for (q = q0; q< qe; q++) *q = '\0';
      Sgets(q0,N-read_so_far,sock);
      // pprint(*lineptr,N);
      for (q = q0; q<qe; q++) if (*q == '\n') break;
      if (q<qe) break;
      read_so_far = N-2;
    }
    N = (N ? 2*N : SGETLINE_INIT_SIZE);
    if (N > SGETLINE_MAX_SIZE) return 0;
    new_line_ptr = (char *) malloc(N);
    assert(new_line_ptr);
    if (*lineptr) {
      memcpy(new_line_ptr,*lineptr,read_so_far);
      free(*lineptr);
    }
    *lineptr = new_line_ptr;
    new_line_ptr = NULL;
  }
  return strlen(*lineptr)+1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
inline double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int talk(Socket *sock,comm_mode &mode) { 
  static char *line = NULL, *buf=NULL;
  static size_t N=0,NN=0;
  int stop;

  if (mode==SEND) {
    printf("> ");
    getline (&line,&N,stdin);
    Sprintf(sock,line);
    if (!strcmp(line,"OVER\n")) mode = RECV;
    stop = !strcmp(line,"STOP\n");
  } else if (mode==RECV) {
    Sgetline(&buf,&NN,sock);
    if (!strcmp(buf,"OVER\n")) mode = SEND;
    printf("-- %s",buf);
    stop = !strcmp(buf,"STOP\n");
  } else assert(0);
  return stop;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Socket *SRVR, *srvr;
int connected=0;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void *wait_connection(void *arg) {
  srvr = Saccept(SRVR);
  connected=1;
  return NULL;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int main(int argc,char **args) {
  comm_mode mode;
  Socket *clnt = NULL;
  pthread_t *thread;
  double inside=0., counter=0.;
  int ierr;
  if(argc>1 && !strcmp(args[1],"-server")) {
    SRVR = Sopen("","s5555");
    assert(SRVR);

    while (1) {
      printf("Waiting connection...\n");
      connected=0;
      ierr =  pthread_create(thread,NULL,&wait_connection,NULL);
      int computed=0;
      while (1) {
	int chunk = 100;
	for (int j=0; j<chunk; j++) {
	  computed++;
	  double x,y;
	  x = drand();
	  y = drand();
	  counter += 1.0;
	  if (x*x+y*y<1.) inside += 1.;
	}
	if (connected) break;
      }
      double cpi = 4.*inside/counter;
      printf("computed %d, current pi %f, error %g\n",
	     computed,cpi,fabs(cpi-M_PI));

      mode = SEND;
      while (!talk(srvr,mode));
      Sclose(srvr);
    }
    
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
