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
  int i,N, *icone_p,node,nread,nnod,nnod2,ndim,ndof;
  double *xnod_p,*data_p;
  Array icone=NULL,xnod=NULL,data=NULL; 
  Group g=NULL;
  char *socket_name_p;
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

  clnt = Sopen("","c5555");
  Sgets(buf,BUFSIZE,clnt);
  sscanf(buf,"nodes %d %d",&ndim,&nnod);
  DXMessage("Got nnod %d, ndim %d",nnod,ndim);
  
  g = DXNewGroup();
  if (!g) goto error;

  xnod = DXNewArray(TYPE_DOUBLE, CATEGORY_REAL, 1,ndim);
  if (!xnod) goto error;
  xnod = DXAddArrayData(xnod, 0, nnod, NULL);
  if (!xnod) goto error;
  xnod_p = (double *)DXGetArrayData(xnod);

  nread = Sreadbytes(clnt,xnod_p,ndim*nnod*sizeof(double));
  if (nread==EOF) goto error;
  g = DXSetMember(g,"nodes",(Object)xnod);
  if (!g) goto error;

  Sgets(buf,BUFSIZE,clnt);
  sscanf(buf,"fields %d %d",&ndof,&nnod2);
  if (nnod!=nnod2) goto error;
  DXMessage("Got ndof %d",ndof);

  data = DXNewArray(TYPE_DOUBLE, CATEGORY_REAL, 1, ndof);
  if (!data) goto error;
  data = DXAddArrayData(data, 0, nnod, NULL);
  if (!data) goto error;
  data_p = (double *)DXGetArrayData(data);
  nread = Sreadbytes(clnt,data_p,ndof*nnod*sizeof(double));
  g = DXSetMember(g,"data",(Object)data);
  if (!g) goto error;

  out[0] = (Object)g;
  Sclose(clnt);

  return OK;

error:
  return ERROR;
}

