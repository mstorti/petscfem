//__INSERT_LICENSE__
//$Id: dxhook.cpp,v 1.9 2003/02/07 23:18:32 mstorti Exp $
#ifdef USE_SSL

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include <src/texthf.h>
#include <src/hook.h>
#include <src/dxhook.h>

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
#define SGETLINE_FACTOR 2
#define SGETLINE_INIT_SIZE 512
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
void dx_hook::init(Mesh &mesh_a,Dofmap &dofmap_a,
			   const char *name_a) {
  int ierr;
  char skthost[10];
  //o TCP/IP port for communicating with DX (5000 < dx_port < 65536). 
  TGETOPTDEF(mesh_a.global_options,int,dx_port,5314);
  TGETOPTDEF_ND(mesh_a.global_options,int,steps,1);
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
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::time_step_pre(double time,int step) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int string2int(string s,int &n) {
  int nread = sscanf(s.c_str(),"%d",&n);
  return (nread!=1);
} 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::
time_step_post(double time,int step,
	       const vector<double> &gather_values) {
  int ierr;
  PetscPrintf(PETSC_COMM_WORLD,
	      "On step %d, step_cntr %d, steps %d\n",
	      step,step_cntr,steps);
  if (steps && step_cntr--) return;
  PetscPrintf(PETSC_COMM_WORLD,"Passed step %d\n",step);

#define BUFSIZE 512
  static char *buf = (char *)malloc(BUFSIZE);
  static size_t Nbuf = BUFSIZE;

  Nodedata *nodedata = mesh->nodedata;
  double *xnod = nodedata->nodedata;
  int ndim = nodedata->ndim;
  int nnod = nodedata->nnod;
  int nu = nodedata->nu;
  PetscPrintf(PETSC_COMM_WORLD,
	      "dx_hook: accepting connections, step %d\n",step);
  if (!MY_RANK) {
    srvr = Saccept(srvr_root);
    assert(srvr);

#if 0
    Sgetline(&buf,&Nbuf,srvr);
#else // DEBUG -----
    Sgets(buf,Nbuf,srvr);
    printf("Got buf %s\n",buf);
#endif
  }
  int Nbuff = Nbuf;
  ierr = MPI_Bcast (&Nbuf, 1, MPI_INT, 0,PETSC_COMM_WORLD);
  if (MY_RANK && Nbuf>Nbuff) {
    free(buf); buf = (char *)malloc(Nbuf);
  }
  ierr = MPI_Bcast (buf, Nbuf, MPI_CHAR, 0,PETSC_COMM_WORLD);
  vector<string> tokens;
  tokenize(buf,tokens);

  // Parse DX options
  int j=0;
  while (1) {
    if (j>=tokens.size()) break;
    if (tokens[j]=="steps") {
      int stepso;
      assert(!string2int(tokens[++j],stepso));
      if (stepso>=0) steps=stepso;
    } else {
      PetscPrintf(PETSC_COMM_WORLD,
		  "Unknown option \"%s\"\n",tokens[j].c_str());
    }
    j++;
  }
  step_cntr = steps-1;
  ierr = MPI_Bcast (&step_cntr, 1, MPI_INT, 0,PETSC_COMM_WORLD);
  
#if 1
  if (!MY_RANK) {
    // Send node coordinates
    Sprintf(srvr,"nodes nodes %d %d\n",ndim,nnod);
    for (int node=0; node<nnod; node++) {
#ifndef DX_USE_FLOATS
      Swrite(srvr,xnod+node*nu,ndim*sizeof(double));
#else
      for (int j=0; j<ndim; j++) {
	float val = (float)*(xnod+node*nu+j);
	Swrite(srvr,&val,sizeof(float));
      }
#endif
    }
  }
  // Send results
  int ndof = dofmap->ndof;
  
  double *state_p = NULL;
  if (!MY_RANK) state_p = new double[ndof*nnod];
  ierr = state2fields(state_p,state(),dofmap,time_data()); assert(!ierr);
  if (!MY_RANK) {
    Sprintf(srvr,"state state %d %d\n",ndof,nnod);
#ifndef DX_USE_FLOATS
    Swrite(srvr,state_p,ndof*nnod*sizeof(double));
#else
    for (int j=0; j<ndof*nnod; j++) {
      float val = (float)*(state_p+j);
      Swrite(srvr,&val,sizeof(float));
    }
#endif
  }
  delete[] state_p;

  // Send connectivities for each elemset
  Darray *elist = mesh->elemsetlist;
  for (int j=0; j<da_length(elist); j++) {
    Elemset *e = *(Elemset **)da_ref(elist,j);
    e->dx(srvr,nodedata,state_p);
  }

  // Send termination signal
  if (!MY_RANK) Sprintf(srvr,"end\n");
#endif

  Sclose(srvr);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::close() {
  Sclose(srvr_root);
}

#endif
