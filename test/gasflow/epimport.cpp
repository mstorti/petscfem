/*
 * Automatically generated on "/tmp/epimport.mb" by DX Module Builder
 */
#include <string>
#include <vector>
#include <map>

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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Converts a line in a list of tokens separated by white space. 
    #tokens# is cleared before the tokenization. 
    @param line (input) line to be tokenized
    @param tokens (output) vector f tokens
*/ 
void tokenize(const char *line,vector<string> &tokens) {
  // Make a local copy (input is read only)
  char *copy = new char[strlen(line)+1];
  // White space pattern
  char spc[] = "[ \t\n]";
  // Clear tokens arg
  tokens.clear();
  // Tokenize using `strtok'
  int j=0;
  while(1) {
    char *token = strtok((j ? NULL : copy),spc);
    tokens.push_back(token);
    j++;
  }
  // clear local copy
  delete[] copy;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define SGETLINE_FACTOR 2
// #define SGETLINE_INIT_SIZE 128
#define SGETLINE_INIT_SIZE 16
#define SGETLINE_MAX_SIZE INT_MAX
/** Reads a line from a socket using the Simple sockets 
    library function #Sgets# but with eventual reallocation, 
    using #malloc#. (This is similar ro the GNU #getline# function). 
    @param lineptr (input/output) the buffer where characters are read. 
    After use, you can free with #free#. 
    @param N_a (input/output) number of bytes initially allocated in #lineptr#
    @param (input) the socket where the line is read. 
    @return number of bytes read */ 
ssize_t Sgetline(char **lineptr, size_t *N_a,Socket *sock) {
  unsigned int &N = *N_a;	// use reference for better readbility
  char * new_line_ptr = NULL, *q, *q0, *qe;
  // At any time the buffer is #N# bytes long and we have
  // read already #read_so_far# bytes.
  int read_so_far=0;
  // Main loop. We read lines with gets until a "\n" is found. If the line
  // has not a "\n" then it should end in "\0\0". 
  while (1) {
    if (N>0) {
      // We read onto lineptr[q0,qe)
      q0 = *lineptr+read_so_far;
      qe = *lineptr+N;
      // Set all to nulls (if a `\n' is left, then we could detect
      // a false line. 
      for (q = q0; q< qe; q++) *q = '\0';
      // Get next part of the line
      Sgets(q0,N-read_so_far,sock);
      // pprint(*lineptr,N);
      // If a newline is found, then we have read the line
      for (q = q0; q<qe; q++) if (*q == '\n') break;
      if (q<qe) break;
      // If not end, then we have read all, except the two nulls at the end
      read_so_far = N-2;
    }
    // Realocate
    N = (N ? 2*N : SGETLINE_INIT_SIZE);
    if (N > SGETLINE_MAX_SIZE) return 0;
    new_line_ptr = (char *) malloc(N);
    assert(new_line_ptr);
    // If already have a buffer, copy to new allocated. 
    if (*lineptr) {
      memcpy(new_line_ptr,*lineptr,read_so_far);
      free(*lineptr);
    }
    // update pointer
    *lineptr = new_line_ptr;
    new_line_ptr = NULL;
  }
  // return number of bytes read
  return strlen(*lineptr)+1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
class GenericError : public string {};
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int string2int(string s,int &n) {
  int nread = sscanf(s.c_str(),"%d",&n);
  return (nread!=1);
} 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
typedef map<string,Array> positions_table_t;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Error build_dx_array(Socket *clnt,int shape,int size, Array &array) {
  array = NULL;
  array = DXNewArray(TYPE_DOUBLE, CATEGORY_REAL, 1,shape);
  if (!array) return ERROR;
  array = DXAddArrayData(array, 0, size, NULL);
  if (!array) return ERROR;
  double *array_p = (double *)DXGetArrayData(array);

  int nread = Sreadbytes(clnt,array_p,shape*size*sizeof(double));
  if (nread==EOF) return ERROR;
  return OK;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
extern "C" Error m_ExtProgImport(Object *in, Object *out) {
  int i,N, *icone_p,node,nread,nnod,nnod2,ndim,ndof,
    nelem,nel;
  double *xnod_p,*data_p;
  Array icone=NULL,xnod=NULL,data=NULL; 
  Group g=NULL;
  char *token;
  String s;
  Type t;
  Socket *clnt = NULL;
  vector<string> tokens;
  string name; 
  positions_table_t positions_table;
  positions_table_t::iterator q, qe=positions_table.end();
  Array array = NULL;
  Error err;
  char spc[] = " \t\n";

#define BUFSIZE 512
  static char *buf = (char *)malloc(BUFSIZE);
  static size_t Nbuf = 0;

  out[0] = NULL;

  char *hostname = DXGetString((String)in[0]); 
  DXMessage("Got hostname: %s",hostname);
  int port;
  if (!in[1]) port = 5314;
  else if (!DXExtractInteger(in[1],&port)) {
    DXSetError(ERROR_DATA_INVALID,
	       "Couldn't find an integer on \"port\" entry");
    goto error;
  }
  
  if (port<=5000 || port>=65536) {
    DXSetError(ERROR_DATA_INVALID,
	       "Invalid port %d, should be in range 5000 < port < 65536",port);
    goto error;
  }
  DXMessage("Got port: %d",port);
  char sktport[20];
  sprintf(sktport,"c%d",port);

  clnt = Sopen(hostname,sktport);
  if (!clnt) {
    DXSetError(ERROR_INTERNAL, "Couldn't open socket");
    return ERROR;
  }
    
  g = DXNewGroup();
  if (!g) goto error;

  while(1) {
    Sgetline(&buf,&Nbuf,clnt);
    DXMessage("Got buf %s",buf);
    tokenize(buf,tokens);
    
    if (tokens[0]=="end") break;
    else if (tokens[0]=="nodes") {
      name = tokens[1];
      if (string2int(tokens[2],ndim)) goto error;
      if (string2int(tokens[3],nnod)) goto error;
      DXMessage("Got nnod %d, ndim %d",nnod,ndim);
      err = build_dx_array(clnt,ndim,nnod,array);
      if(err!=OK) return err;
      positions_table[name] = array;
    };
  }

  for (q=positions_table.begin(); q!=qe; q++) {
    g = DXSetMember(g,(char *)(q->first.c_str()),
		    (Object)q->second);
  }

#if 0
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

#endif
  out[0] = (Object)g;

  if (!clnt) Sclose(clnt);
  clnt = NULL;
  return OK;

error:
  Sclose(clnt);
  delete[] hostname;
  hostname = NULL;
  return ERROR;
} 
