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
  int N=10, *icone_p, j,k,base, elem=0;
  Array icone=NULL; 
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
  /* Set `connections' component */
  f = DXSetComponentValue(f,"connections",(Object)icone); if (!f) goto error;
  f = DXEndField(f); if (!f) goto error;
  out[0] = (Object)f;

  return OK;

error:
  return ERROR;
}

