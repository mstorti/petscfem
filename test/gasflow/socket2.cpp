//__INSERT_LICENSE__
// $Id: socket2.cpp,v 1.3 2003/02/03 15:32:07 mstorti Exp $
#define _GNU_SOURCE
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <HDR/sockets.h>

enum comm_mode { SEND, RECV };

#define BUFSIZE 200

void chomp(char *s) { s[strlen(s)-1] = '\0'; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int talk(Socket *sock,comm_mode &mode) { 
  static char *line = NULL;
  static size_t N=0;
  char buf[100], *buff;
  int stop;

  if (mode==SEND) {
    printf("> ");
    getline (&line,&N,stdin);
    Sprintf(sock,line);
    if (!strcmp(line,"OVER\n")) mode = RECV;
    stop = !strcmp(line,"STOP\n");
  } else if (mode==RECV) {
    buff = Sgets(buf,BUFSIZE,sock);
    assert(buff);
    if (!strcmp(buf,"OVER\n")) mode = SEND;
    printf("-- \"%s\"",buf);
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
