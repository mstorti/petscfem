//__INSERT_LICENSE__
// $Id: epimport.cpp,v 1.6 2003/02/16 01:42:01 mstorti Exp $
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
#define USE_SSL
#include <src/util3.h>

static Error traverse(Object *, Object *);
static Error doLeaf(Object *, Object *);

#define DXassert(cond) 								\
if(!(cond)) {									\
  DXMessage("Assertion \"%s\" failed at %s:%d",#cond,__FILE__,__LINE__);	\
  goto error;									\
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
  int rank,size,nnod;
  vector<int> shape;
  Array array;
  Object dx_object() { return (Object)array; }
  State(int k,vector<int> &s,int z,int d,Array a) 
    : rank(k), shape(s), size(z), nnod(d), array(a) {}
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
#define DX_USE_FLOATS

#ifdef DX_USE_FLOATS
#define DX_SCALAR_TYPE TYPE_FLOAT
#define DX_SCALAR_TYPE_DECL float
#else
#define DX_SCALAR_TYPE TYPE_DOUBLE
#define DX_SCALAR_TYPE_DECL double
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Error build_dx_array(Socket *clnt,int shape,int length, Array &array) {
  array = NULL;
  array = DXNewArray(DX_SCALAR_TYPE, CATEGORY_REAL, 1,shape);
  if (!array) return ERROR;
  array = DXAddArrayData(array, 0, length, NULL);
  if (!array) return ERROR;
  DX_SCALAR_TYPE_DECL *array_p = (DX_SCALAR_TYPE_DECL *)DXGetArrayData(array);

  int nread = Sreadbytes(clnt,array_p,shape*length*sizeof(DX_SCALAR_TYPE_DECL));
  if (nread==EOF) return ERROR;
  return OK;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Error build_dx_array_v(Socket *clnt,int rank,
		       vector<int> shape, int size, 
		       int length, Array &array) {
  array = NULL;
  array = DXNewArrayV(DX_SCALAR_TYPE, CATEGORY_REAL, rank, shape.begin());
  if (!array) return ERROR;
  array = DXAddArrayData(array, 0, length, NULL);
  if (!array) return ERROR;
  DX_SCALAR_TYPE_DECL *array_p = (DX_SCALAR_TYPE_DECL *)DXGetArrayData(array);

  int nread = Sreadbytes(clnt,array_p,size*length*sizeof(DX_SCALAR_TYPE_DECL));
  if (nread==EOF) return ERROR;
  return OK;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Error build_dx_array_int(Socket *clnt,int shape,int length, Array &array) {
  array = NULL;
  array = DXNewArray(TYPE_INT, CATEGORY_REAL, 1,shape);
  if (!array) return ERROR;
  array = DXAddArrayData(array, 0, length, NULL);
  if (!array) return ERROR;
  int *array_p = (int *)DXGetArrayData(array);

  int nread = Sreadbytes(clnt,array_p,shape*length*sizeof(int));
  if (nread==EOF) return ERROR;
  return OK;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class AutoString {
private:
  char *s;
  int n;
public:
  AutoString() : s(NULL) {}
  ~AutoString() { free(s); }
  char *str() const { return s; }
  int *N() { return &n; }
  void resize(int m) { 
    if (m>n) {
      char *new_s = (char *)malloc(m);
      if (n>0) {
	strcpy(new_s,s);
	free(s);
      }
      n=m;
    }
  }
  void clear() {
    if (n>0) { n=0; free(s); }
  }
  void cat(AutoString &s) { 
    // Cats s at the rear of this string
    // not coded yet
  }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
extern "C" Error m_ExtProgImport(Object *in, Object *out) {
  int i,N, *icone_p,node,nread,nnod,nnod2,ndim,ndof,
    nelem, nel, cookie, fields_n, arrays_n;

  double *xnod_p,*data_p;
  Array icone=NULL,xnod=NULL,data=NULL,step_comp_o=NULL; 
  Group alist=NULL,flist=NULL;
  char *token;
  String s;
  Type t;
  Socket *clnt = NULL;
  vector<string> tokens;
  string name; 
  DXObjectsTable dx_objects_table;
  DXObjectsTable::iterator q,r,qe;
  Array array = NULL;
  Error ierr;
  char spc[] = " \t\n";

#define BUFSIZE 512
  static char *buf = (char *)malloc(BUFSIZE);
  static size_t Nbuf = BUFSIZE;

  out[0] = NULL;
  out[1] = NULL;
  out[2] = NULL;

  // Inputs
  int in_index = 0;
  Object steps_o          = in[in_index++];
  Object hostname_o       = in[in_index++];
  Object port_o           = in[in_index++];
  Object options_o        = in[in_index++];
  Object step_o           = in[in_index++];
  Object state_file_o     = in[in_index++];

  int steps, port, step, step_comp;
  char *options, *hostname, *state_file;
  DXMessage("before processing steps");
  if (!steps_o) steps = -1;
  else if (!DXExtractInteger(steps_o,&steps)) {
    DXSetError(ERROR_DATA_INVALID,
	       "Couldn't find an integer on \"port\" entry");
    goto error;
  }
  if (steps<-1) steps=-1;

  options = DXGetString((String)options_o); 
  if (!options) options = "";

  hostname = DXGetString((String)hostname_o); 
  if (!hostname) hostname = "localhost";

  if (!port_o) port = 5314;
  else if (!DXExtractInteger(port_o,&port)) {
    DXSetError(ERROR_DATA_INVALID,
	       "Couldn't find an integer on \"port\" entry");
    goto error;
  }
  
  if (port<=5000 || port>=65536) {
    DXSetError(ERROR_DATA_INVALID,
	       "Invalid port %d, should be in range 5000 < port < 65536",port);
    goto error;
  }

  if (!step_o) step = -1;
  else if (!DXExtractInteger(step_o,&step)) {
    DXSetError(ERROR_DATA_INVALID,
	       "Couldn't find an integer on \"step\" entry");
    goto error;
  }

  state_file = DXGetString((String)state_file_o); 
  if (!state_file) state_file="<no-state>";

  DXMessage("Got steps %d, hostname \"%s\", port %d, "
	    "step %d, state_file \"%s\", other options \"%s\"",
	    steps,hostname,port,step,state_file,options);


  char sktport[20];
  sprintf(sktport,"c%d",port);

  clnt = Sopen(hostname,sktport);
  if (!clnt) {
    DXSetError(ERROR_INTERNAL, "Couldn't open socket");
    return ERROR;
  }

  DXMessage("Sending steps %d %s",steps,options);
  Sprintf(clnt,"steps %d step %d state_file %s %s\n",
	  steps,step,state_file,options);

  Sgetline(&buf,&Nbuf,clnt);
  tokenize(buf,tokens);
  if (tokens[0]!="step") {
    DXMessage("Not \"step\" found in first line, got %s",buf);
    goto error;
  }
  if (string2int(tokens[1],step_comp)) goto error;
  DXMessage("Got computation step %d",step_comp);

  while(1) {
    Sgetline(&buf,&Nbuf,clnt);
    DXMessage("Got line \"%s\"",buf);
    tokenize(buf,tokens);

    if (tokens[0]=="end") break;
    //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    else if (tokens[0]=="nodes") {
      string &name = tokens[1];
      if (string2int(tokens[2],ndim)) goto error;
      if (string2int(tokens[3],nnod)) goto error;
      if (string2int(tokens[4],cookie)) goto error;
      ierr = build_dx_array(clnt,ndim,nnod,array);
      if(ierr!=OK) return ierr;
      ierr = dx_objects_table.load_new(name,new Nodes(ndim,nnod,array));
      if(ierr!=OK) return ierr;
      DXMessage("Got new \"Nodes\" name %s, ptr %p, ndim %d, nnod %d",
		name.c_str(),array,ndim,nnod);
      Sprintf(clnt,"nodes_OK %d\n",cookie);
      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    } else if (tokens[0]=="state") {
      DXMessage("Got state line %s",buf);
      int rank,dim,tk=1,size=1;
      vector<int> shape;
      string &name = tokens[tk++];
      if (string2int(tokens[tk++],rank)) goto error;
      for (int j=0; j<rank; j++) {
	if (string2int(tokens[tk++],dim)) goto error;
	size *= dim;
	shape.push_back(dim);
      }
      if (string2int(tokens[tk++],nnod)) goto error;
      if (string2int(tokens[tk++],cookie)) goto error;

      ierr = build_dx_array_v(clnt,rank,shape,size,nnod,array);
      if(ierr!=OK) return ierr;
      ierr = dx_objects_table.
	load_new(name,new State(rank,shape,size,nnod,array));
      if(ierr!=OK) return ierr;
      DXMessage("Got new \"State\" name %s, ptr %p, rank %d, dims (",
		name.c_str(),array,rank);
      for (int j=0; j<rank; j++) DXMessage(" %d",shape[j]);
      DXMessage("), size %d, nnod %d",size,nnod);
      Sprintf(clnt,"state_OK %d\n",cookie);
      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    } else if (tokens[0]=="elemset") {
      DXMessage("Got line %s",buf);
      string &name = tokens[1];
      string &dx_type = tokens[2];
      if (string2int(tokens[3],nel)) goto error;
      if (string2int(tokens[4],nelem)) goto error;
      if (string2int(tokens[5],cookie)) goto error;
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
      Sprintf(clnt,"elemset_OK %d\n",cookie);
      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    } else if (tokens[0]=="field") {
      // Get components 
      Object positions,connections,data;
      string &name = tokens[1];
      string &pname = tokens[2];	// name of "positions" component
      ierr = dx_objects_table.get_positions(pname,positions);
      if(ierr!=OK) return ierr;
      string &cname = tokens[3];
      ierr = dx_objects_table.get_connections(cname,connections);
      if(ierr!=OK) return ierr;
      string &fname = tokens[4];
      ierr = dx_objects_table.get_state(fname,data);
      if(ierr!=OK) return ierr;
      if (string2int(tokens[5],cookie)) goto error;
      Sprintf(clnt,"field_OK %d\n",cookie);

      // Build new field
      Field field = DXNewField();
      if (!field) goto error;
      field = DXSetComponentValue(field,"positions",(Object)positions); 
      if (!field) goto error;
      field = DXSetComponentValue(field,"connections",(Object)connections); 
      if (!field) goto error;
      field = DXSetComponentValue(field,"data",(Object)data); 
      if (!field) goto error;

      field = DXEndField(field); if (!field) goto error;

      // Load new field in table
      ierr = dx_objects_table.load_new(name,new DXField(pname,cname,fname,field));
      if(ierr!=OK) return ierr;

      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    } else if (tokens[0]=="fields_auto") {

      if (string2int(tokens[1],cookie)) goto error;
      DXMessage("Got \"fields_auto\" directive, cookie %d",cookie);

      Object positions,connections,data;
      string p("nodes");
      ierr = dx_objects_table.get_positions(p,positions);
      if(ierr!=OK) return ierr;

      qe=dx_objects_table.end();
      for (q=dx_objects_table.begin(); q!=qe; q++) {
	Elemset *elemset = dynamic_cast<Elemset *>(q->second);
	if (!elemset) continue;
	Object connections = elemset->dx_object();
	for (r=dx_objects_table.begin(); r!=qe; r++) {
	  State *state = dynamic_cast<State *>(r->second);
	  if (!state) continue;
	  Object data = state->dx_object();

	  Field field = DXNewField();
	  if (!field) goto error;
	  field = DXSetComponentValue(field,"positions",(Object)positions); 
	  if (!field) goto error;
	  field = DXSetComponentValue(field,"connections",connections); 
	  if (!field) goto error;
	  field = DXSetComponentValue(field,"data",data); 
	  if (!field) goto error;

	  field = DXEndField(field); if (!field) goto error;

	  // Load new field in table
	  string n("nodes");
	  string cname(q->first);
	  string dname(r->first);
	  DXField *dxf = new DXField(n,cname,dname,field);
	  string fname = cname + "_" + dname;
	  DXMessage("Creating field %s",fname.c_str());
	  ierr = dx_objects_table.load_new(fname,dxf);
	  if(ierr!=OK) return ierr;
	}
      }      

      Sprintf(clnt,"fields_auto_OK %d\n",cookie);

    //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    } else {
      DXSetError(ERROR_INTERNAL,
		 "Received bad line %s",tokens[0].c_str());
      goto error;
    }
  }

  fields_n=0; arrays_n=0;
  DXMessage("Before building alist and flist");
  qe = dx_objects_table.end();
  for (q=dx_objects_table.begin(); q!=qe; q++) {
    DXField *field = dynamic_cast<DXField *>(q->second);
    Object o;
    char *name = (char *)(q->first.c_str());
    if (field) {
      if (!fields_n++) {
	flist = DXNewGroup();
	if (!flist) goto error;
      }
      o = (Object)field->dx_object();
      flist = DXSetMember(flist,name,o);
      if (!flist) goto error;
      DXMessage("Field %d, ptr. %p, member name \"%s\"",fields_n,o,name);
    } else {
      if (!arrays_n++) {
	alist = DXNewGroup();
	if (!alist) goto error;
      }
      o = (Object)q->second->dx_object();
      alist = DXSetMember(alist,name,o);
      if (!alist) goto error;
      DXMessage("Array %d, ptr. %p, member name \"%s\"",arrays_n,o,name);
    }
  }

  step_comp_o = DXMakeInteger(step_comp);

  if(alist) out[0] = (Object)alist;
  if(flist) out[1] = (Object)flist;
  out[2] = (Object)step_comp_o;

  if (!clnt) Sclose(clnt);
  clnt = NULL;
  return OK;

error:
  Sclose(clnt);
  delete[] hostname;
  hostname = NULL;
  return ERROR;
} 
