#define _GNU_SOURCE
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int make_socket (uint16_t port) {
  int sock;
  struct sockaddr_in name;

  /* Create the socket. */
  sock = socket (PF_INET, SOCK_STREAM, 0);
  if (sock < 0)
    {
      perror ("socket");
      exit (EXIT_FAILURE);
    }

  /* Give the socket a name. */
  name.sin_family = AF_INET;
  name.sin_port = htons (port);
  name.sin_addr.s_addr = htonl (INADDR_ANY);
  if (bind (sock, (struct sockaddr *) &name, sizeof (name)) < 0) {
      perror ("bind");
      exit (EXIT_FAILURE);
    }
  return sock;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void init_sockaddr (struct sockaddr_in *name,
               const char *hostname,
               uint16_t port) {
  struct hostent *hostinfo;

  name->sin_family = AF_INET;
  name->sin_port = htons (port);
#if 1
  hostinfo = gethostbyname (hostname);
  if (hostinfo == NULL) {
      fprintf (stderr, "Unknown host %s.\n", hostname);
      exit (EXIT_FAILURE);
  }
  name->sin_addr = *(struct in_addr *) hostinfo->h_addr;
#else
  name->sin_addr.s_addr = INADDR_LOOPBACK;
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void write_to_server (int filedes) {
  int nbytes;
  char *message = "This is the message...";
  nbytes = write (filedes, message, strlen (message) + 1);
  if (nbytes < 0) {
    perror ("write");
    exit (EXIT_FAILURE);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int read_from_client (int filedes) {
#define MAXMSG 200
  char buffer[MAXMSG];
  int nbytes;
  
  nbytes = read (filedes, buffer, MAXMSG);
  if (nbytes < 0) {
    /* Read error. */
    perror ("read");
    exit (EXIT_FAILURE);
  } else if (nbytes == 0) return -1; // End-of-file.
  else {
    // Data read.
    fprintf (stderr, "Server: got message: `%s'\n", buffer);
    return 0;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int read_data(int s,char *buf,int n) {
  /*s =  connected socket */
  /* buf = pointer to the buffer */
  /* n = number of characters (bytes) we want */
  int bcount, br;/* counts bytes read, bytes read this pass */
  
  bcount= 0;
  br= 0;
  while (bcount < n) {             /* loop until full buffer */
    if ((br= read(s,buf,n-bcount)) > 0) {
      bcount += br;                /* increment byte counter */
      buf += br;                   /* move buffer ptr for next read */
    }
    if (br < 0) return(-1);        /* signal an error to the caller */
  }
  return(bcount);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int write_data(int s,const char *buf,int n) {
  /* s = connected socket */
  /* buf = pointer to the buffer */
  /* n = number of characters (bytes) we want */
  int bcount, br;/* counts bytes read, bytes read this pass */
  
  bcount= 0;
  br= 0;
  while (bcount < n) {             /* loop until full buffer */
    if ((br= write(s,buf,n-bcount)) > 0) {
      bcount += br;                /* increment byte counter */
      buf += br;                   /* move buffer ptr for next read */
    }
    if (br < 0) return(-1);        /* signal an error to the caller */
  }
  return(bcount);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void chomp(char *s) { s[strlen(s)-1] = '\0'; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define PORT 5555
#define SERVERHOST "spider"
#define BUFSIZE 100

enum comm_mode { SEND, RECV };

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int talk(int sock,comm_mode &mode) { 
  static char *line = NULL;
  static size_t N=0;
  char buf2[100];
  if (mode==SEND) {
    printf("> ");
    getline (&line,&N,stdin);
    assert(strlen(line)<BUFSIZE);
    chomp(line);
    strcpy(buf2,line);
    // printf("sending \"%s\"",buf2);
    write_data(sock,buf2,BUFSIZE);
    if (!strcmp(line,"OVER")) mode = RECV;
  } else if (mode==RECV) {
    read_data(sock,buf2,BUFSIZE);
    if (!strcmp(buf2,"OVER")) mode = SEND;
    printf("-- %s\n",buf2);
  } else assert(0);
  return !strcmp(buf2,"STOP");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int main(int argc,char **args) {
  int s, sock;
  // fd_set active_fd_set, read_fd_set;
  sockaddr_in servername, clientname;
  comm_mode mode;
  // char buf[BUFSIZE];
  const char *buf = "Hello socket!\n\0";
  if(argc>1 && !strcmp(args[1],"-server")) {
    sock = make_socket (PORT);
    printf("server: trace 0\n");
    if (listen (sock, 1) < 0) {
	perror ("listen");
	exit (EXIT_FAILURE);
    }
    printf("server: trace 0.1\n");
    size_t size = sizeof(clientname);
    s = accept(sock,(struct sockaddr *)&clientname,&size);
    if (s < 0) {
      perror ("accept");
      exit (EXIT_FAILURE);
    }
    mode = SEND;
    while (!talk(s,mode));
    close(s);
    exit(EXIT_SUCCESS);

  } else {
    sock = socket (PF_INET, SOCK_STREAM, 0);
    if (sock < 0) {
      perror("socket (client)");
      exit(EXIT_FAILURE);
    }
    /* Connect to the server. */
    init_sockaddr(&servername, SERVERHOST, PORT);
    if (0 > connect (sock, (struct sockaddr *) &servername, sizeof (servername))) {
      perror ("connect (client)");
      exit (EXIT_FAILURE);
    }
    printf("client: trace 0\n");
    mode = RECV;
    while (!talk(sock,mode));
    close(sock);
    exit(EXIT_SUCCESS);
  }

}
