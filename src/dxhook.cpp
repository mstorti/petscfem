//__INSERT_LICENSE__
//$Id: dxhook.cpp,v 1.1 2003/02/03 15:52:56 mstorti Exp $

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include <src/texthf.h>
#include <src/hook.h>
#include <src/dxhook.h>

#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>

#define PF_DX_PORT 5555

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
void dx_hook::init(Mesh &mesh,Dofmap &dofmap,
			   const char *name_a) {
  int s, sock;
  sockaddr_in servername, clientname;
  const char *buf = "Hello socket!\n\0";

  options = new TextHashTableFilter(mesh.global_options);
  options->push("dx");

  sock = make_socket (PF_DX_PORT);
  if (listen (sock, 1) < 0) {
    perror ("listen");
    exit (EXIT_FAILURE);
  }
  // printf("server: trace 0.1\n");
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
}
