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
    int, int *);

#if defined (__cplusplus) || defined (c_plusplus)
extern "C"
#endif
Error
m_ExtProgImport(Object *in, Object *out)
{
  int i;

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


  /*
   * Output "value" is a value; set up an appropriate array
   */
  out[0] = (Object)DXNewArray(TYPE_INT, CATEGORY_REAL, 0, 0);
  if (! out[0])
    goto error;
  if (! DXAddArrayData((Array)out[0], 0, 1, NULL))
    goto error;
  memset(DXGetArrayData((Array)out[0]), 0, 1*DXGetItemSize((Array)out[0]));


  /*
   * Call the hierarchical object traversal routine
   */
  if (!traverse(in, out))
    goto error;

  return OK;

error:
  /*
   * On error, any successfully-created outputs are deleted.
   */
  for (i = 0; i < 1; i++)
  {
    if (in[i] != out[i])
      DXDelete(out[i]);
    out[i] = NULL;
  }
  return ERROR;
}


static Error
traverse(Object *in, Object *out)
{
  switch(DXGetObjectClass(in[0]))
  {
    case CLASS_FIELD:
    case CLASS_ARRAY:
    case CLASS_STRING:
      /*
       * If we have made it to the leaf level, call the leaf handler.
       */
      if (! doLeaf(in, out))
  	     return ERROR;

      return OK;

    case CLASS_GROUP:
    {
      int   i, j;
      int   memknt;
      Class groupClass  = DXGetGroupClass((Group)in[0]);

      DXGetMemberCount((Group)in[0], &memknt);


       /*
        * Create new in and out lists for each child
        * of the first input. 
        */
        for (i = 0; i < memknt; i++)
        {
          Object new_in[1], new_out[1];

         /*
          * For all inputs that are Values, pass them to 
          * child object list.  For all that are Field/Group, get 
          * the appropriate decendent and place it into the
          * child input object list.
          */

          /* input "socket_name" is Value */
          new_in[0] = in[0];

         /*
          * For all outputs that are Values, pass them to 
          * child object list.  For all that are Field/Group,  get
          * the appropriate decendent and place it into the
          * child output object list.  Note that none should
          * be NULL (unlike inputs, which can default).
          */

          /* output "value" is Value */
          new_out[0] = out[0];

          if (! traverse(new_in, new_out))
            return ERROR;

        }
      return OK;
    }

    case CLASS_XFORM:
    {
      int    i, j;
      Object new_in[1], new_out[1];


      /*
       * Create new in and out lists for the decendent of the
       * first input.  For inputs and outputs that are Values
       * copy them into the new in and out lists.  Otherwise
       * get the corresponding decendents.
       */

      /* input "socket_name" is Value */
      new_in[0] = in[0];

      /*
       * For all outputs that are Values, copy them to 
       * child object list.  For all that are Field/Group,  get
       * the appropriate decendent and place it into the
       * child output object list.  Note that none should
       * be NULL (unlike inputs, which can default).
       */

      /* output "value" is Value */
      new_out[0] = out[0];

      if (! traverse(new_in, new_out))
        return ERROR;

      return OK;
    }

    case CLASS_SCREEN:
    {
      int    i, j;
      Object new_in[1], new_out[1];


      /*
       * Create new in and out lists for the decendent of the
       * first input.  For inputs and outputs that are Values
       * copy them into the new in and out lists.  Otherwise
       * get the corresponding decendents.
       */

      /* input "socket_name" is Value */
      new_in[0] = in[0];


      /*
       * For all outputs that are Values, copy them to 
       * child object list.  For all that are Field/Group,  get
       * the appropriate decendent and place it into the
       * child output object list.  Note that none should
       * be NULL (unlike inputs, which can default).
       */

       /* output "value" is Value */
       new_out[0] = out[0];

       if (! traverse(new_in, new_out))
         return ERROR;

       return OK;
     }

     case CLASS_CLIPPED:
     {
       int    i, j;
       Object new_in[1], new_out[1];


       /* input "socket_name" is Value */
       new_in[0] = in[0];


      /*
       * For all outputs that are Values, copy them to 
       * child object list.  For all that are Field/Group,  get
       * the appropriate decendent and place it into the
       * child output object list.  Note that none should
       * be NULL (unlike inputs, which can default).
       */

       /* output "value" is Value */
       new_out[0] = out[0];

       if (! traverse(new_in, new_out))
         return ERROR;

       return OK;
     }

     default:
     {
       DXSetError(ERROR_BAD_CLASS, "encountered in object traversal");
       return ERROR;
     }
  }
}

static int
doLeaf(Object *in, Object *out)
{
  int i, result=0;
  Array array;
  Field field;
  Pointer *in_data[1], *out_data[1];
  int in_knt[1], out_knt[1];
  Type type;
  Category category;
  int rank, shape;
  Object attr, src_dependency_attr = NULL;
  char *src_dependency = NULL;
  int p_knt = -1;
  int c_knt = -1;

  if (DXGetObjectClass(in[0]) == CLASS_FIELD)
  {

    field = (Field)in[0];

    if (DXEmptyField(field))
      return OK;

    /* 
     * Determine the dependency of the source object's data
     * component.
     */
    src_dependency_attr = DXGetComponentAttribute(field, "data", "dep");
    if (! src_dependency_attr)
    {
      DXSetError(ERROR_MISSING_DATA, "\"socket_name\" data component is missing a dependency attribute");
      goto error;
    }

    if (DXGetObjectClass(src_dependency_attr) != CLASS_STRING)
    {
      DXSetError(ERROR_BAD_CLASS, "\"socket_name\" dependency attribute");
      goto error;
    }

    src_dependency = DXGetString((String)src_dependency_attr);

    array = (Array)DXGetComponentValue(field, "positions");
    if (! array)
    {
      DXSetError(ERROR_BAD_CLASS, "\"socket_name\" contains no positions component");
      goto error;
    }

    DXGetArrayInfo(array, &p_knt, NULL, NULL, NULL, NULL);
    /* 
     * If there are connections, get their count so that
     * connections-dependent result arrays can be sized.
     */
    array = (Array)DXGetComponentValue(field, "connections");
    if (array)
        DXGetArrayInfo(array, &c_knt, NULL, NULL, NULL, NULL);
  }
  /*
   * If the input argument is not NULL then we get the 
   * data array: either the object itself, if its an 
   * array, or the data component if the argument is a field
   */
  if (! in[0])
  {
    array = NULL;
    in_data[0] = NULL;
    in_knt[0] = NULL;
  }
  else
  {
    if (DXGetObjectClass(in[0]) == CLASS_ARRAY)
    {
      array = (Array)in[0];
    }
    else if (DXGetObjectClass(in[0]) == CLASS_STRING)
    {
      in_data[0] = (Pointer *)DXGetString((String)in[0]);
      in_knt[0] = 1;
    }
    else
    {
      if (DXGetObjectClass(in[0]) != CLASS_FIELD)
      {
        DXSetError(ERROR_BAD_CLASS, "\"socket_name\" should be a field");
        goto error;
      }

      array = (Array)DXGetComponentValue((Field)in[0], "data");
      if (! array)
      {
        DXSetError(ERROR_MISSING_DATA, "\"socket_name\" has no data component");
        goto error;
      }

      if (DXGetObjectClass((Object)array) != CLASS_ARRAY)
      {
        DXSetError(ERROR_BAD_CLASS, "data component of \"socket_name\" should be an array");
        goto error;
      }
    }


    if (DXGetObjectClass(in[0]) != CLASS_STRING)    {
       DXGetArrayInfo(array, &in_knt[0], &type, &category, &rank, &shape);
       if (type != TYPE_STRING || category != CATEGORY_REAL ||
             !((rank == 0) || ((rank == 1)&&(shape == 1))))
       {
         DXSetError(ERROR_DATA_INVALID, "input \"socket_name\"");
         goto error;
       }

       in_data[0] = DXGetArrayData(array);
       if (! in_data[0])
          goto error;

    }
  }
  if (! out[0])
  {
    DXSetError(ERROR_INTERNAL, "Value output 0 (\"value\") is NULL");
    goto error;
  }
  if (DXGetObjectClass(out[0]) != CLASS_ARRAY)
  {
    DXSetError(ERROR_INTERNAL, "Value output 0 (\"value\") is not an array");
    goto error;
  }

  array = (Array)out[0];

  DXGetArrayInfo(array, &out_knt[0], &type, &category, &rank, &shape);
  if (type != TYPE_INT || category != CATEGORY_REAL ||
      !((rank == 0) || ((rank == 1)&&(shape == 1))))
  {
    DXSetError(ERROR_DATA_INVALID, "Value output \"value\" has bad type");
    goto error;
  }

  /*
   * Now get the pointer to the contents of the array
   */
  out_data[0] = DXGetArrayData(array);
  if (! out_data[0])
    goto error;

  /*
   * Call the user's routine.  Check the return code.
   */
  result = ExtProgImport_worker(
    in_knt[0], (char *)in_data[0],
    out_knt[0], (int *)out_data[0]);

  if (! result)
     if (DXGetError()==ERROR_NONE)
        DXSetError(ERROR_INTERNAL, "error return from user routine");

  /*
   * In either event, clean up allocated memory
   */

error:
  return result;
}

int
ExtProgImport_worker(
    int socket_name_knt, char *socket_name_data,
    int value_knt, int *value_data)
{
  /*
   * The arguments to this routine are:
   *
   *
   * The following are inputs and therefore are read-only.  The default
   * values are given and should be used if the knt is 0.
   *
   * socket_name_knt, socket_name_data:  count and pointer for input "socket_name"
   *                   no default value given.
   *
   *  The output data buffer(s) are writable.
   *  The output buffer(s) are preallocated based on
   *     the dependency (positions or connections),
   *     the size of the corresponding positions or
   *     connections component in the first input, and
   *     the data type.
   *
   * value_knt, value_data:  count and pointer for output "value"
   */

  /*
   * User's code goes here
   */
  *value_data = 3141;
     
  /*
   * successful completion
   */
   return 1;
     
  /*
   * unsuccessful completion
   */
error:
   return 0;
  
}
