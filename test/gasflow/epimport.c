/*
 * Automatically generated on "/tmp/epimport.mb" by DX Module Builder
 */

/* define your pre-dx.h include file for inclusion here*/ 
#ifdef PRE_DX_H
#include "ExtProgImport_predx.h"
#endif
#include "dx/dx.h"
/* define your post-dx.h include file for inclusion here*/ 
#ifdef POST_DX_H
#include "ExtProgImport_postdx.h"
#endif
#include <HDR/sockets.h>

static Error traverse(Object *, Object *);
static Error doLeaf(Object *, Object *);

/*
 * Declare the interface routine.
 */
int
ExtProgImport_worker(
    int, char *,
    int, double *);

#if defined (__cplusplus) || defined (c_plusplus)
extern "C"
#endif
Error
m_ExtProgImport(Object *in, Object *out)
{
  int i,N, *icone_p, j,k,base, elem=0, node,nread,
    nnod,ndim;
  float *xnod_p,x,y,*data_p,r2,c;
  Array icone=NULL,xnod=NULL,data=NULL; 
  char *socket_name_p;
  Field f=NULL; 
  String s;
  Type t;
  Socket *clnt;
#define BUFSIZE 512
  char buf[BUFSIZE];

  /*
   * Initialize all outputs to NULL
   */
  out[0] = NULL;

  /*
   * Error checks: required inputs are verified.
   */

  /* Parameter "socket_name" is required. */
  if (in[0] == NULL) {
    DXSetError(ERROR_MISSING_DATA, "\"socket_name\" must be specified");
    return ERROR;
  }

  DXMessage("trace 0");
  clnt = Sopen("","c5555");
  DXMessage("trace 1");
  Sgets(buf,BUFSIZE,clnt);
  DXMessage("got %s",buf);
  DXMessage("trace 2");
  sscanf(buf,"nodes %d %d",&ndim,&nnod);
  DXMessage("Got nnod %d, ndim %d",nnod,ndim);
  Sclose(clnt);
  DXMessage("trace 3");

  return OK;

error:
  return ERROR;
}

