//__INSERT_LICENSE__
//$Id: inviscid.cpp,v 1.3 2003/01/01 15:11:15 mstorti Exp $
#define _GNU_SOURCE

extern int MY_RANK,SIZE;

#include "./inviscid.h"
#include "./fifo.h"

extern int MY_RANK,SIZE;
#define NDIM 2    
#define NDOF (NDIM+1)    

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class ext_node {
public:
  double x[NDIM], u[NDIM];
};

vector<int> visc_nodes,ext_nodes;
typedef map<int,ext_node> ext_map;
ext_map coup_vals;

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
    
    FILE *fid = fopen("ext.coupling_nodes.tmp","r");
    int nread;
    while (1) {
      int vn,en;
      nread = fscanf(fid,"%d %d",&vn,&en);
      if (nread!=2) break;
      visc_nodes.push_back(vn);
      ext_nodes.push_back(en);
      coup_vals[en] = ext_node();      
    }
    assert(nread==EOF);

    fid = fopen("cylin.nod.tmp","r");
    double x[NDIM];
    int node=0;
    ext_map::iterator q, qe=coup_vals.end();
    while (1) {
      nread = fscanf(fid,"%lf %lf",&x[0],&x[1]);
      if(nread!=NDIM) break;
      node++;
      q = coup_vals.find(node);
      if (q!=qe) {
	q->second.x[0] = x[0];
	q->second.x[1] = x[1];
      }
    }
    assert(nread==EOF);
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
public:
  coupling() { }
  void init(TextHashTable *thash);
  double eval(double);
  ~coupling() { }
};

DEFINE_EXTENDED_AMPLITUDE_FUNCTION(coupling);
