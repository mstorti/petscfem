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
/* #include <math.h> */

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
  int i;
  int N=40, *icone_p, j,k,base, elem=0, node;
  float *xnod_p,x,y,*data_p,r2;
  Array icone=NULL,xnod=NULL,data=NULL; 
  Field f=NULL; 

  /*
   * Initialize all outputs to NULL
   */
  out[0] = NULL;

  /*
   * Error checks: required inputs are verified.
   */

  /* Parameter "socket_name" is required. */
  if (in[0] == NULL)
  {
    DXSetError(ERROR_MISSING_DATA, "\"socket_name\" must be specified");
    return ERROR;
  }

  f = DXNewField();
  if (!f) goto error;

  /* ====================================================== */
  icone = DXNewArray(TYPE_INT, CATEGORY_REAL, 1,4);
  if (!icone) goto error;
  icone = DXAddArrayData(icone, 0, N*N, NULL);
  if (!icone) goto error;
  icone_p = (int *)DXGetArrayData(icone);
  /* Define connectivities    */
  for (j=0; j<N; j++) {
    for (k=0; k<N; k++) {
      base = j*(N+1)+k;
      *icone_p++ = base;
      *icone_p++ = base+1;
      *icone_p++ = base+N+1;
      *icone_p++ = base+N+2;
    }
  }
  icone = (Array)DXSetStringAttribute((Object)icone,
				      "element type","quads"); if (!icone) goto error;
  /* Set `connections' component */
  f = DXSetComponentValue(f,"connections",(Object)icone); if (!f) goto error;
  /* attribute "element type" string "quads" */

  /* ====================================================== */
  xnod = DXNewArray(TYPE_FLOAT, CATEGORY_REAL, 1,2);
  if (!xnod) goto error;
  xnod = DXAddArrayData(xnod, 0, (N+1)*(N+1), NULL);
  if (!xnod) goto error;
  xnod_p = (float *)DXGetArrayData(xnod);

  data = DXNewArray(TYPE_FLOAT, CATEGORY_REAL, 0);
  if (!data) goto error;
  data = DXAddArrayData(data, 0, (N+1)*(N+1), NULL);
  if (!data) goto error;
  data_p = (float *)DXGetArrayData(data);

  /* Define positions and results */
  for (j=0; j<=N; j++) {
    x = ((float) j)/((float) N);
    for (k=0; k<=N; k++) {
      y = ((float) k)/((float) N);
      *xnod_p++ = x;
      *xnod_p++ = y;
#define SQ(a) ((a)*(a))
      r2 = SQ(x-0.5) + SQ(y-0.5);
      *data_p++ = 0.1/(0.1+r2);
    }
  }
  /* Set `connections' component */
  f = DXSetComponentValue(f,"positions",(Object)xnod); if (!f) goto error;
  f = DXSetComponentValue(f,"data",(Object)data); if (!f) goto error;

  /* ====================================================== */
  f = DXEndField(f); if (!f) goto error;
  out[0] = (Object)f;

  return OK;

error:
  return ERROR;
}

