//__INSERT_LICENSE__
// $Id: epimport.cpp,v 1.8 2003/02/07 13:19:26 mstorti Exp $
#include <string>
#include <vector>
#include <map>
#include <strstream>

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
  strcpy(copy,line);
  // White space pattern
  char spc[] = "[ \t\n]";
  // Clear tokens arg
  tokens.clear();
  // Tokenize using `strtok'
  int j=0;
  while(1) {
    char *token = strtok((j ? NULL : copy),spc);
    if (!token) break;
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
int string2int(string s,int &n) {
  int nread = sscanf(s.c_str(),"%d",&n);
  return (nread!=1);
} 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class DXObject {
public:
  virtual Object dx_object() { return NULL; }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class Nodes : public DXObject {
public:
  int ndim,nnod;
  Array array;
  Object dx_object() { return (Object)array; }
  Nodes(int m,int d,Array a) : ndim(m), nnod(d), array(a) {}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class State : public DXObject {
public:
  int ndof,nnod;
  Array array;
  Object dx_object() { return (Object)array; }
  State(int f,int d,Array a) : ndof(f), nnod(d), array(a) {}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class Elemset : public DXObject {
public:
  int nel,nelem;
  string dx_type;
  Array array;
  Object dx_object() { return (Object)array; }
  Elemset(int l,int m,string &dxt,Array a) : nel(l), nelem(m), 
    dx_type(dxt), array(a) {}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class DXField : public DXObject {
public:
  string positions,connections,data;
  Field field;
  Object dx_object() { return (Object)field; }
  DXField(string &p, string &c, string &d,Field f) : positions(p),
  connections(c), data(d), field(f) {}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class DXObjectsTable : public map<string,DXObject *> {
private:
  static int field_sfx_max;
public:
  ~DXObjectsTable() {
    map<string,DXObject *>::iterator q,qe;
    for (q=begin(); q!=end(); q++) {
      delete q->second;
      q->second = NULL;
    }
  }
  int get_positions(string &name,Object &object);
  int get_connections(string &name,Object &object);
  int get_state(string &name,Object &object);
  int load_new(string &name,DXObject *dx_object);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int DXObjectsTable::field_sfx_max = 1000;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int DXObjectsTable::load_new(string &name,DXObject *dx_object) {
  static char *buf = NULL;
  static int Nbuf=0;
  if (!dx_object) {
    DXSetError(ERROR_INTERNAL,
	       "Attempts to load null object");
    return ERROR;
  }
  iterator q,qe;
  string new_name;
  if (find(name)!=end()) {
    int j;
    for (j=0; j < field_sfx_max; j++) {
      Nbuf = asprintf(&buf,"%s_%d",name.c_str(),j);
      assert(Nbuf>=0);
      new_name = string(buf);
      if (find(new_name)!=end()) break;
    }
    if (j<field_sfx_max) 
      DXMessage("renaming field \"%s\" -> \"%s\" to avoid collision",
		name.c_str(), new_name.c_str());
    else {
      DXSetError(ERROR_INTERNAL, "Can't rename field \"%s\"");
      return ERROR;
    }
  } else new_name = name;
  (*this)[new_name] = dx_object;
  return OK;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
int DXObjectsTable::get_positions(string &name,Object &object) {
  iterator q = find(name);
  if (q==end()) {
    DXSetError(ERROR_DATA_INVALID,
	       "Can't find object \"%s\"",name.c_str());
    return ERROR;
  }
  Nodes *nodes = dynamic_cast<Nodes *>(q->second);
  if (!nodes) {
    DXSetError(ERROR_DATA_INVALID,
	       "Can't convert object \"%s\" to type Nodes",name.c_str());
    return ERROR;
  }
  object = nodes->dx_object();
  return OK;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
int DXObjectsTable::get_connections(string &name,Object &object) {
  iterator q = find(name);
  if (q==end()) {
    DXSetError(ERROR_DATA_INVALID,
	       "Can't find object \"%s\"",name.c_str());
    return ERROR;
  }
  Elemset *elemset = dynamic_cast<Elemset *>(q->second);
  if (!elemset) {
    DXSetError(ERROR_DATA_INVALID,
	       "Can't convert object \"%s\" to type Elemset",name.c_str());
    return ERROR;
  }
  object = elemset->dx_object();
  return OK;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
int DXObjectsTable::get_state(string &name,Object &object) {
  iterator q = find(name);
  if (q==end()) {
    DXSetError(ERROR_DATA_INVALID,
	       "Can't find object \"%s\"",name.c_str());
    return ERROR;
  }
  State *state = dynamic_cast<State *>(q->second);
  if (!state) {
    DXSetError(ERROR_DATA_INVALID,
	       "Can't convert object \"%s\" to type State",name.c_str());
    return ERROR;
  }
  object = state->dx_object();
  return OK;
}

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
Error build_dx_array_int(Socket *clnt,int shape,int size, Array &array) {
  array = NULL;
  array = DXNewArray(TYPE_INT, CATEGORY_REAL, 1,shape);
  if (!array) return ERROR;
  array = DXAddArrayData(array, 0, size, NULL);
  if (!array) return ERROR;
  int *array_p = (int *)DXGetArrayData(array);

  int nread = Sreadbytes(clnt,array_p,shape*size*sizeof(int));
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
  DXObjectsTable dx_objects_table;
  DXObjectsTable::iterator q, qe=dx_objects_table.end();
  Array array = NULL;
  Error ierr;
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

  Array p,c,d;
  Field f;
  while(1) {
    Sgetline(&buf,&Nbuf,clnt);
    DXMessage("Got buf %s",buf);
    tokenize(buf,tokens);
    DXMessage("After tokenize");
    
    if (tokens[0]=="end") break;
    //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    else if (tokens[0]=="nodes") {
      string &name = tokens[1];
      if (string2int(tokens[2],ndim)) goto error;
      if (string2int(tokens[3],nnod)) goto error;
      ierr = build_dx_array(clnt,ndim,nnod,array);
      if(ierr!=OK) return ierr;
      ierr = dx_objects_table.load_new(name,new Nodes(ndim,nnod,array));
      if(ierr!=OK) return ierr;
      DXMessage("Got new \"Nodes\" name %s, ptr %p, ndim %d, nnod %d",
		name.c_str(),array,ndim,nnod);
      p = array;
    //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    } else if (tokens[0]=="state") {
      name = tokens[1];
      if (string2int(tokens[2],ndof)) goto error;
      if (string2int(tokens[3],nnod)) goto error;
      ierr = build_dx_array(clnt,ndof,nnod,array);
      if(ierr!=OK) return ierr;
      ierr = dx_objects_table.load_new(name,new State(ndim,nnod,array));
      if(ierr!=OK) return ierr;
      DXMessage("Got new \"State\" name %s, ptr %p, ndof %d, nnod %d",
		name.c_str(),array,ndof,nnod);
      d = array;
    //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    } else if (tokens[0]=="elemset") {
      string &name = tokens[1];
      string &dx_type = tokens[2];
      if (string2int(tokens[3],nel)) goto error;
      if (string2int(tokens[4],nelem)) goto error;
      ierr = build_dx_array_int(clnt,nel,nelem,array);
      if(ierr!=OK) return ierr;
      array = (Array)
	DXSetStringAttribute((Object)array,
			     "element type",(char *)dx_type.c_str()); 
      if (!array) goto error;
      ierr = dx_objects_table
	.load_new(name,new Elemset(nel,nelem,dx_type,array));
      if(ierr!=OK) return ierr;
      DXMessage("Got new \"Elemset\" name %s, ptr %p, nel %d, nelem %d",
		name.c_str(),array,nel,nelem);
      c = array;
    //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    } else if (tokens[0]=="field") {
#if 0
      // Get components 
      Object positions,connections,data;
      string &name = tokens[1];
      string &pname = tokens[2];	// name of "positions" component
      ierr = dx_objects_table.get_positions(pname,positions);
      if(ierr!=OK) return ierr;
      DXMessage("trace 0");
      string &cname = tokens[3];
      ierr = dx_objects_table.get_connections(cname,connections);
      DXMessage("trace 1");
      if(ierr!=OK) return ierr;
      string &fname = tokens[4];
      ierr = dx_objects_table.get_state(fname,data);
      DXMessage("trace 2");
      if(ierr!=OK) return ierr;

      // Build new field
      Field field = DXNewField();
      if (!field) goto error;
      field = DXSetComponentValue(field,"positions",(Object)positions); 
      DXMessage("trace 3, positions %p",positions);
      if (!field) goto error;
      field = DXSetComponentValue(field,"connections",(Object)connections); 
      if (!field) goto error;
      DXMessage("trace 4, connections %p",connections);
      field = DXSetComponentValue(field,"data",(Object)data); 
      if (!field) goto error;
      DXMessage("trace 5, data %p",data);
      field = DXEndField(field); if (!field) goto error;
      DXMessage("trace 6");

      // Load new field in table
      ierr = dx_objects_table.load_new(name,new DXField(pname,cname,fname,field));
      if(ierr!=OK) return ierr;
      DXMessage("trace 7");
#endif
      DXMessage("trace 8");
      Field field = DXNewField();
      if (!field) goto error;
      DXMessage("trace 9");
      field = DXSetComponentValue(field,"positions",(Object)p); 
      if (!field) goto error;
      DXMessage("trace 10");
      field = DXSetComponentValue(field,"connections",(Object)c); 
      if (!field) goto error;
      DXMessage("trace 11");
      field = DXSetComponentValue(field,"data",(Object)d); 
      if (!field) goto error;
      DXMessage("trace 12");
      field = DXEndField(field); if (!field) goto error;
      f = field;
      DXMessage("trace 13");

    //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    } else {
      DXSetError(ERROR_INTERNAL,
		 "Received bad line %s",tokens[0].c_str());
      goto error;
    }
  }

#if 0
  g = DXNewGroup();
  if (!g) goto error;

  for (q=dx_objects_table.begin(); q!=qe; q++) {
    g = DXSetMember(g,(char *)(q->first.c_str()),
		    (Object)q->second->dx_object());
  }

  out[0] = (Object)g;
#endif
  out[0] = (Object)f;

  if (!clnt) Sclose(clnt);
  clnt = NULL;
  return OK;

error:
  Sclose(clnt);
  delete[] hostname;
  hostname = NULL;
  return ERROR;
} 
