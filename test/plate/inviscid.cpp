//__INSERT_LICENSE__
//$Id: inviscid.cpp,v 1.2 2002/12/30 03:06:31 mstorti Exp $
#define _GNU_SOURCE

extern int MY_RANK,SIZE;

#include "./inviscid.h"
#include "./fifo.h"

extern int MY_RANK,SIZE;


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void coupling_visc_hook::init(Mesh &mesh,Dofmap &dofmap,
		     TextHashTableFilter *options,const char *name) {
  if (!MY_RANK) {
    printf("VISCOUS: Opening fifos for communicating with INVISCID.\n");

    visc2inv = fopen("visc2inv.fifo","w");
    assert(visc2inv);
    setvbuf(visc2inv,NULL,_IOLBF,0);

    inv2visc = fopen("inv2visc.fifo","r");
    assert(inv2visc);

    printf("VISCOUS: Done.\n");
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void coupling_visc_hook::time_step_pre(double t,int step) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void coupling_visc_hook::time_step_post(double time,int step,
			       const vector<double> &gather_values) {
  fprintf(visc2inv,"computed_step %d\n",step,time);
}

DL_GENERIC_HOOK(coupling_visc_hook);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void coupling_inv_hook::init(Mesh &mesh,Dofmap &dofmap,
		     TextHashTableFilter *options,const char *name) {
  assert(SIZE==1); // Not implemented in parallel yet
  if (!MY_RANK) {
    printf("INVISCID: Opening fifos for communicating with VISCOUS.\n");

    inv2visc = fopen("inv2visc.fifo","w");
    assert(inv2visc);
    setvbuf(inv2visc,NULL,_IOLBF,0);

    visc2inv = fopen("visc2inv.fifo","r");
    assert(visc2inv);

    printf("INV2VISC: Done.\n");
    
#define NDIM 2    
    FILE *fid = fopen("cylin.nod.tmp","r");
    double x[NDIM];
    int nread;
    while (1) {
      nread = fscanf(fid,"%lf %lf",&x[0],&x[1]);
      if(nread!=NDIM) break;
    }
    assert(nread==EOF);
    fclose(fid);
    nnod = xnod.size()/NDIM;

#define NDOF (NDIM+1)    
    u.resize(nnod*NDOF);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void coupling_inv_hook::time_step_pre(double t,int step) { 
  double step_sent = read_doubles(visc2inv,"computed_step");
  assert(int(step_sent)==step);
  FILE *fid = fopen("cylin.state.tmp","r");
  int nread;
  for (int j=0; j<nnod*NDOF; j++) {
    nread = fscanf(fid,"%lf",&u[j]);
    if(nread!=1) break;
  }
  assert(nread==EOF);
  fclose(fid);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void coupling_inv_hook::time_step_post(double time,int step,
			       const vector<double> &gather_values) {
  fprintf(inv2visc,"computed_step %d\n",step);
}

DL_GENERIC_HOOK(coupling_inv_hook);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class coupling {
private:
  int *nodes;
  double *vals;
public:
  coupling() : nodes(NULL), vals(NULL) {}
  void init(TextHashTable *thash);
  double eval(double);
  ~coupling() {
    delete[] nodes; nodes=NULL;
    delete[] vals; vals=NULL;
  }
};

DEFINE_EXTENDED_AMPLITUDE_FUNCTION(coupling);
