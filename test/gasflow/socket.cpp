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
#if 0
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
#define PORT 5555
#define SERVERHOST "spider"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int main(int argc,char **args) {
  int sock;
  fd_set active_fd_set, read_fd_set;
  sockaddr_in servername, clientname;
  if(argc>1 && !strcmp(args[1],"-server")) {
    sock = make_socket (PORT);
    printf("server: trace 0\n");
    if (listen (sock, 1) < 0) {
	perror ("listen");
	exit (EXIT_FAILURE);
    }
    printf("server: trace 0.1\n");
    FD_ZERO (&active_fd_set);
    FD_SET (sock, &active_fd_set);
    read_fd_set = active_fd_set;
    if (select (FD_SETSIZE, &read_fd_set, NULL, NULL, NULL) < 0) {
      perror ("select");
      exit (EXIT_FAILURE);
    }
    if (FD_ISSET (sock, &read_fd_set)) {
      /* Connection request on original socket. */
      size_t size = sizeof (clientname);
      if (accept (sock, (struct sockaddr *) &clientname, &size) < 0) {
	perror ("accept");
	exit (EXIT_FAILURE);
      }
      printf("server: trace 1\n");
      fprintf (stderr,
	       "Server: connect from host %s, port %hd.\n",
	       inet_ntoa (clientname.sin_addr),
	       ntohs (clientname.sin_port));
      printf("server: trace 2\n");
    }
    /* Initialize the set of active sockets. */
    // FD_ZERO (&active_fd_set);
    // FD_SET (sock, &active_fd_set);
  } else{
    sock = socket (PF_INET, SOCK_STREAM, 0);
    if (sock < 0) {
      perror ("socket (client)");
      exit (EXIT_FAILURE);
    }
    /* Connect to the server. */
    init_sockaddr(&servername, SERVERHOST, PORT);
    printf("client: trace 0\n");
    write_to_server (sock);
    printf("client: trace 1\n");
    close (sock);
    printf("client: trace 2\n");
    exit (EXIT_SUCCESS);
  }
}
