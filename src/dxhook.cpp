//__INSERT_LICENSE__
//$Id: dxhook.cpp,v 1.41 2003/06/08 13:10:43 mstorti Exp $

#include <src/debug.h>
#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include <src/util3.h>
#include <src/texthf.h>
#include <src/hook.h>
#include <src/dxhook.h>
#include <src/autostr.h>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <src/sockbuff.h>

#ifdef USE_SSL

#include <HDR/sockets.h>

extern int MY_RANK, SIZE;

// I didn't found a definition for this
#ifndef IPPORT_MAX
#define IPPORT_MAX 65536
#endif

// It seems that DX hangs when using doubles for coordinates and
// data in the moment of making `DXEndField()' (it seems that internally
// it happens in th moment of doing `DXBoundingBox()'. So that
// I will use only floats. 
#define DX_USE_FLOATS

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FieldGenList::~FieldGenList() {
  for (int j=0; j<size(); j++) {
    delete (*this)[j];
    (*this)[j] = NULL;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class FieldGenDefault : public FieldGen {
  int ndof;
  string name;
public:
  void init(int ndof_a,TextHashTable* options,char *name_a) { 
    ndof=ndof_a; 
    name = string(name_a);
  }
  int n() { return 1; }
  void field(int j,string &name,vector<int> &shape) {
    name = "state";
    shape.clear();
    shape.push_back(ndof);
  }
  void values(int j,vector<double> &in,vector<double> &out) {
    for (int k=0; k<ndof; k++) out[k] = in[k];
  }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class FieldGenLine : public FieldGen {
  int ndof;
  string name;
  struct entry {
    string name;
    vector<int> rank;
    vector<int> dof_list;
    int size;
  };
  vector<entry> field_list;
public:
  void parse(const char *line);
  void init(int ndof_a,TextHashTable* options,char *name_a) { 
    ndof=ndof_a; 
    name = string(name_a);
  }
  int n() { return field_list.size(); }
  void field(int j,string &name,vector<int> &rank) {
    const entry &e = field_list[j];
    name = e.name;
    rank = e.rank;
  }
  void values(int jf,vector<double> &in,vector<double> &out) {
    const vector<int> &dof_list = field_list[jf].dof_list;
    for (int j=0; j<dof_list.size(); j++) out[j] = in[dof_list[j]];
  }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FieldGenLine::parse(const char *line) {
  vector<string> tokens;
  tokenize(line,tokens);
  int ierr;
  int jtok=0;
  int ntoks = tokens.size();
  while (jtok < ntoks) {
    entry e;
    int rank;
    ierr = string2int(tokens[jtok++],rank);
    assert(!ierr);
    assert(rank>=0);
    assert(jtok+rank <= ntoks);
    e.size = 1;
    for (int j=0; j<rank; j++) {
      int dim;
      ierr = string2int(tokens[jtok++],dim);
      assert(!ierr);
      assert(dim>0);
      e.rank.push_back(dim);
      e.size *= dim;
    }
    assert(jtok + e.size <= ntoks);
    e.dof_list.resize(e.size);
    for (int j=0; j<e.size; j++) {
      int dof;
      ierr = string2int(tokens[jtok++],dof);
      assert(!ierr);
      assert(dof>=0 && dof<ndof);
      e.dof_list[j] = dof;
    }
    assert(jtok<ntoks);
    e.name = tokens[jtok++];
    field_list.push_back(e);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
dx_hook::dx_hook() : 
#ifdef USE_PTHREADS
  connection_state_m(not_launched),
  connection_state_master(not_launched) , 
#endif
  options(NULL), srvr_root(NULL), step_cntr(0), steps(0) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::init(Mesh &mesh_a,Dofmap &dofmap_a,
			   const char *name_a) {
  int ierr;
  char skthost[10];
  //o TCP/IP port for communicating with DX ({\tt 5000 < dx\_port < 65536}). 
  TGETOPTDEF(mesh_a.global_options,int,dx_port,5314);
  //o Initial value for the {\tt steps} parameter. 
  TGETOPTDEF(mesh_a.global_options,int,dx_steps,1);
  steps = dx_steps;
  if (!MY_RANK) {
    PETSCFEM_ASSERT(dx_port>IPPORT_USERRESERVED && 
		    dx_port<IPPORT_MAX,"\"dx_port\" number must be in the range\n"
		    "IPPORT_USERRESERVED < dx_port < IPPORT_MAX, dx_port = %d\n"
		    "[current values are: IPPORT_USERRESERVED = %d\n, IPPORT_MAX  = %d]\n",
		    dx_port, IPPORT_USERRESERVED, IPPORT_MAX);
    printf("dx_hook: starting socket at port: %d\n",dx_port);
    sprintf(skthost,"S%d",dx_port);
    srvr_root = Sopen("",skthost);
    assert(srvr_root);
    PetscPrintf(PETSC_COMM_WORLD,"Done.\n");
  }
  mesh = &mesh_a;
  dofmap = &dofmap_a;
  
  TextHashTable *go = mesh_a.global_options;
  int dx_split_state_flag=0;

  //o Read states from file instead of computing them . Normally
  //  this is done to analyze a previous run. 
  TGETOPTDEF(go,int,dx_read_state_from_file,0);

  //o Generates DX fields by combination of the input fields
  TGETOPTDEF_S(go,string,dx_split_state,);
  if (dx_split_state!="") {
    dx_split_state_flag = 1;
    FieldGenLine *fgl = new FieldGenLine;
    fgl->init(dofmap->ndof,mesh_a.global_options,"split_state");
    fgl->parse(dx_split_state.c_str());
    field_gen_list.push_back(fgl);
  }

  //o Generates a DX field with the whole state (all ndof fields)
  TGETOPTDEF(go,int,dx_state_all_fields,!dx_split_state_flag);
  if (dx_state_all_fields) {
    FieldGenDefault *fg = new FieldGenDefault;
    fg->init(dofmap->ndof,mesh_a.global_options,"all_fields");
    field_gen_list.push_back(fg);
  }

  //o Auto generate states by combining elemsets with fields
  TGETOPTDEF_ND(go,int,dx_auto_combine,0);

  if (dx_read_state_from_file) {
    if (steps==0) steps = 1;
    int step=0;
    while(1) {
      send_state(step++,&dx_hook::build_state_from_file);
      if (record==-1) {
	PetscFinalize();
	exit(0);
      }
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::time_step_pre(double time,int step) {}

#define PF_DBG(name,format)					\
PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] " 		\
                        #name "  " #format "\n",MY_RANK,name);	\
PetscSynchronizedFlush(PETSC_COMM_WORLD); 
#define PF_DBG_INT(name)  PF_DBG(name,%d) 
#define PF_DBG_DBL(name)  PF_DBG(name,%f) 

#ifdef USE_PTHREADS
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
static void *wait_connection(void *arg) {
  dx_hook *hook = (dx_hook *)arg;
  return hook->wait_connection();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void *dx_hook::wait_connection() {
  srvr = Saccept(srvr_root);
  connection_state_master = connected;
  return NULL;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
dx_hook::connection_state_t 
dx_hook::connection_state() {
  if (!MY_RANK) connection_state_m = connection_state_master;
  ierr = MPI_Bcast (&connection_state_m, 1, MPI_INT, 0,PETSC_COMM_WORLD);
  return connection_state_m;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::set_connection_state(connection_state_t s) {
  connection_state_m = connection_state_master = s;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::re_launch_connection() {
  set_connection_state(not_connected);
  if (!MY_RANK) {
    ierr = pthread_create(&thread,NULL,&::wait_connection,this);
    assert(!ierr);
  }
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int dx_hook::build_state_from_state(double *state_p) {
  return state2fields(state_p,state(),dofmap,time_data()); 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class GenericError : public string { 
public:
  GenericError(char *s) : string(s) { }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int dx_hook::build_state_from_file(double *state_p) {
#if 1 // debug:=
  int record_local = (record>=0 ? record : 0);
#define record record_local
#endif
  PetscPrintf(PETSC_COMM_WORLD,
	      "dx_hook: reading state from file %s, record %d\n",
	      state_file.c_str(),record);
  if (!MY_RANK) {
    FILE *fid = fopen(state_file.c_str(),"r");
    if(!fid) throw GenericError("Can't open file.");
    int base = record*nnod*ndof;
    for (int j=0; j<(record+1)*nnod*ndof; j++) {
      double val;
      int nread = fscanf(fid,"%lf",&val);
      if(nread!=1) {
	fclose(fid);
	throw GenericError("Can't read line.");
      }
      if (j>base) state_p[j-base] = val;
    }
    fclose(fid);
  }
  return 0;
#undef record
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::send_state(int step,build_state_fun_t build_state_fun) try {
#define sock srvr
  int cookie, cookie2, dx_step;
  if (steps && step_cntr--) return;
  if (!steps) {
#ifdef USE_PTHREADS
    if(connection_state() == not_launched) {
      re_launch_connection();
      sleep(1);
    }
    if (connection_state()==not_connected) {
      PetscPrintf(PETSC_COMM_WORLD,
		  "dx_hook: step %d, detached mode, no connection found\n",step);
      return;
    } else if (connection_state()==connected) {
      void *retval;
      if (!MY_RANK) {
	ierr = pthread_join(thread,&retval);
	assert(!ierr);
      }
      set_connection_state(not_launched);
    } else assert(0);
#else
    PetscPrintf(PETSC_COMM_WORLD,
		"Can't have asynchronous communication (steps=0) with DX \n"
		"because this version is not compiled with threads.\n"
		"Setting steps=1\n");
    steps=1;
#endif
  } else {
    if (!MY_RANK) {
      PetscPrintf(PETSC_COMM_WORLD,
		  "dx_hook: accepting connections, step %d\n",step);
      srvr = Saccept(srvr_root);
      assert(srvr);
    }
  }
#define BUFSIZE 512
  static char *buf = (char *)malloc(BUFSIZE);
  static size_t Nbuf = BUFSIZE;

  vector<string> tokens;
  Nodedata *nodedata = mesh->nodedata;
  double *xnod = nodedata->nodedata;
  ndim = nodedata->ndim;
  nnod = nodedata->nnod;
  ndof = dofmap->ndof;
  int nu = nodedata->nu;
  int stepso;

  // Process DX options. 
  if (!MY_RANK) {
    Sgetline(&buf,&Nbuf,srvr);
    tokenize(buf,tokens);

    // Parse DX options
    int j=0;
    while (1) {
      if (j>=tokens.size()) break;
      if (tokens[j]=="steps") {
	assert(!string2int(tokens[++j],stepso));
      } else if (tokens[j]=="step") {
	assert(!string2int(tokens[++j],dx_step));
      } else if (tokens[j]=="state_file") {
	state_file = tokens[++j];
      } else if (tokens[j]=="record") {
	assert(!string2int(tokens[++j],record));
      } else {
	printf("Unknown option \"%s\"\n",tokens[j].c_str());
      }
      j++;
    }
    printf("dx_hook: Got steps %d, dx_step %d, state_file %s, record %d\n",
	   steps,dx_step,state_file.c_str(),record);
  }

  // Options are read in master and
  // each option is sent to the slaves with MPI_Bcast
  ierr = MPI_Bcast (&stepso, 1, MPI_INT, 0,PETSC_COMM_WORLD);
  ierr = MPI_Bcast (&dx_step, 1, MPI_INT, 0,PETSC_COMM_WORLD);
  ierr = string_bcast(state_file,0,PETSC_COMM_WORLD);
  ierr = MPI_Bcast (&record, 1, MPI_INT, 0,PETSC_COMM_WORLD);

  if (record==-1) throw GenericError("Received record=-1, stop.");

  if (stepso>=0 && stepso!=steps) {
    PetscPrintf(PETSC_COMM_WORLD,
		"dx_hook: changed \"steps\" %d -> %d from DX\n",
		steps,stepso);
    if (stepso==0) set_connection_state(not_launched);
    steps=stepso;
  }
  
  step_cntr = steps-1;
  
  SocketBuffer<float> sbuff(srvr);
  if (!MY_RANK) {
    // Send node coordinates
    cookie = rand();
    Sprintf(srvr,"step %d\n",step);
    Sprintf(srvr,"nodes nodes %d %d %d\n",ndim,nnod,cookie);
    // printf("sending nodes nodes %d %d %d\n",ndim,nnod,cookie);
    for (int node=0; node<nnod; node++)
      for (int j=0; j<ndim; j++) sbuff.put((float)*(xnod+node*nu+j));
    sbuff.flush();
    CHECK_COOKIE(nodes);
  }

  // Send results
  double *state_p=NULL;
  dvector<double> state_v(ndof*nnod);
  if (!MY_RANK) {
    state_v.resize(ndof*nnod);
    assert(state_v.chunks_n()==1);
    state_p = state_v.buff();
  }
  if ((this->*build_state_fun)(state_p)) throw GenericError("Can't build state");
  
  if (!MY_RANK) {
    cookie = rand();
    FieldGenList::iterator qp, qe = field_gen_list.end();
    for (qp=field_gen_list.begin(); qp!=qe; qp++) {
      FieldGen *q = *qp;
      int nf = q->n();
      for (int jf=0; jf<nf; jf++) {
	vector<int> shape;
	string name;
	q->field(jf,name,shape);
	int rank=shape.size(), size=1;
	AutoString buff;
	// Sends name and rank/shape of entity
	buff.sprintf("state %s %d",name.c_str(),rank);

	// Sends the list of dimensions and total size
	for (int jd=0; jd<rank; jd++) {
	  buff.cat_sprintf(" %d",shape[jd]);
	  size *= shape[jd];
	}
	
	// Send number of nodes and cookie
	buff.cat_sprintf(" %d %d",nnod,cookie);
	// printf("sending state line \"%s\"\n",buff.str());
	Sprintf(srvr,"%s\n",buff.str());

	// Send values 
	vector<double> in(ndof),out(size);
	for (int j=0; j<nnod; j++) {
	  double *base_node = state_p+j*ndof;
	  for (int l=0; l<ndof; l++) in[l] = *(base_node+l);
	  q->values(jf,in,out);
	  for (int l=0; l<size; l++) sbuff.put((float)out[l]);
	}
	sbuff.flush();
	CHECK_COOKIE(state);
      }
    }
  }

  // Send connectivities for each elemset
  Darray *elist = mesh->elemsetlist;
  for (int j=0; j<da_length(elist); j++) {
    Elemset *e = *(Elemset **)da_ref(elist,j);
    e->dx(srvr,nodedata,state_p);
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Send signal for automatically pairing connections and fields
  if (dx_auto_combine && !MY_RANK) {
    cookie = rand();
    Sprintf(srvr,"fields_auto %d\n",cookie);
    CHECK_COOKIE(fields_auto);
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Send termination signal
  if (!MY_RANK) {
    Sprintf(srvr,"end\n");
    Sclose(srvr);
    srvr = NULL;
  }
#ifdef USE_PTHREADS
  if(!steps && connection_state() == not_launched) re_launch_connection();
#endif
} catch(GenericError e) {
  if (!MY_RANK && srvr) {
    Sprintf(srvr,"end\n");
    Sclose(srvr);
    printf("%s\n",e.c_str());
  }
  PetscFinalize();
  exit(0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::
time_step_post(double time,int step,
	       const vector<double> &gather_values) {
  send_state(step,&dx_hook::build_state_from_state);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::close() {
#ifdef USE_PTHREADS
  // I'm not too much sure how to do this correctly
  if (connection_state()==connected) {
    // Here we should shut down the connection
    // sending some message to the client
    if (!MY_RANK) {
      Sprintf(srvr,"end\n");
      Sclose(srvr);
    }
    set_connection_state(not_connected);
  }
  if (connection_state()==not_connected) {
    pthread_cancel(thread);
    set_connection_state(not_launched);
  }
  assert(connection_state()==not_launched);
#endif
  if (!MY_RANK) Sclose(srvr_root); 
}

#endif
