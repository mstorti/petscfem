/*
 * Automatically generated on "/tmp/epimport.mb" by DX Module Builder
 */
#include <string>

// `or' and `string' are used in DX so that one possibility is to
// use a namespace for DX or to remap this colliding names with
// `#defines'
#define or __or__
#define string __string__
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
#undef or
#undef string

static Error traverse(Object *, Object *);
static Error doLeaf(Object *, Object *);

#define DXassert(cond) 								\
if(!(cond)) {									\
  DXMessage("Assertion \"%s\" failed at %s:%d",#cond,__FILE__,__LINE__);	\
  goto error;									\
}

#if defined (__cplusplus) || defined (c_plusplus)
extern "C"
#endif
Error
m_ExtProgImport(Object *in, Object *out) {
  int i,N, *icone_p,node,nread,nnod,nnod2,ndim,ndof,
    nelem,nel;
  double *xnod_p,*data_p;
  Array icone=NULL,xnod=NULL,data=NULL; 
  Group g=NULL;
  char *socket_name_p, *token;
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
  socket_name_p = DXGetString((String)in[0]); 
  char *hostname = new char[strlen(socket_name_p)+1];
  if (!socket_name_p) goto error;
  DXMessage("got socket_name_p: %s",socket_name_p);
  int port;
  nread = sscanf(socket_name_p,"%s:%d",hostname,port);
  DXMessage("got nread: %d",nread);
  if (nread!=2) {
    DXSetError(ERROR_DATA_INVALID,
	       "Couldn't parse hostname:port");
    goto error;
  }
  if (!hostname) goto error;
  if (port<=5000 || port>=65536) {
    DXSetError(ERROR_DATA_INVALID,
	       "Invalid port %d, should be in range 5000 < port < 65536",port);
    goto error;
  }
  DXMessage("Got hostname: %s, port: %d",hostname,port);
  char sktport[20];
  sprintf(sktport,"c%d",port);

  clnt = Sopen(hostname,sktport);
  if (!clnt) {
    DXSetError(ERROR_INTERNAL, "Couldn't open socket");
    return ERROR;
  }
    
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

  while(1) {
    char spc[] = " \t\n";
    Sgets(buf,BUFSIZE,clnt);

    token = strtok(buf,spc);
    if(!strcmp(token,"end")) break;
    DXassert(!strcmp(token,"icone"));

    token = strtok(NULL,spc);
    nread = sscanf(token,"%d",&nelem);
    DXassert(nread==1);

    token = strtok(NULL,spc);
    nread = sscanf(token,"%d",&nel);
    DXassert(nread==1);

    token = strtok(NULL,spc);
    DXassert(token);
    DXassert(strlen(token)>0);
    string ename(token);

    token = strtok(NULL,spc);
    DXassert(token);
    DXassert(strlen(token)>0);
    string etype(token);

    DXMessage("Got elemset \"%s\", nelem %d, nel %d, type \"%s\"\n",
	      ename.c_str(),nelem,nel,etype.c_str());

    icone = DXNewArray(TYPE_INT, CATEGORY_REAL, 1, nel);
    if (!icone) goto error;
    icone = DXAddArrayData(icone, 0, nelem, NULL);
    if (!icone) goto error;
    icone_p = (int *)DXGetArrayData(icone);
    nread = Sreadbytes(clnt,icone_p,nelem*nel*sizeof(int));
    g = DXSetMember(g,(char *)ename.c_str(),(Object)icone);
    if (!g) goto error;
  }

  out[0] = (Object)g;
  Sclose(clnt);

  return OK;

error:
  delete[] hostname;
hostname = NULL;
return ERROR;
}

