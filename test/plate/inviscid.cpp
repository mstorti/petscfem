//__INSERT_LICENSE__
//$Id: inviscid.cpp,v 1.5 2003/01/01 16:16:36 mstorti Exp $
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
  // viscous node, external potential node
  int vn,pn;
  // position, velocity, potential at node
  double x[NDIM], u[NDIM], phi;
};

typedef map<int,int> ext_map;
ext_map coup_vals;
vector<ext_node> ext_node_data;

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

    visc2inv = fopen("visc2inv.fifo","r");
    assert(visc2inv);

    inv2visc = fopen("inv2visc.fifo","w");
    assert(inv2visc);
    setvbuf(inv2visc,NULL,_IOLBF,0);

    printf("INV2VISC: Done.\n");
    
    FILE *fid = fopen("ext.coupling_nodes.tmp","r");
    int nread, ext_node_indx=0;
    while (1) {
      int vn,pn;
      nread = fscanf(fid,"%d %d",&vn,&pn);
      if (nread!=2) break;
      ext_node_data.push_back(ext_node());
      ext_node_data[ext_node_indx].vn = vn;
      ext_node_data[ext_node_indx].pn = pn;
      coup_vals[pn] = ext_node_indx;      
      ext_node_indx++;
    }
    assert(nread==EOF);
    fclose(fid);

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
	ext_node_indx = q->second;
	ext_node_data[ext_node_indx].x[0] = x[0];
	ext_node_data[ext_node_indx].x[1] = x[1];
      }
    }
    assert(nread==EOF);
    fclose(fid);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void coupling_inv_hook::time_step_pre(double t,int step) { 
  double step_sent = read_doubles(visc2inv,"computed_step");
  assert(int(step_sent)==step);
  printf("INVISCID: received computed_step flag OK, step %d\n",step);
  FILE *fid = fopen("cylin.state.tmp","r");
  double u[NDIM];
  int node=0;
  ext_map::iterator q, qe=coup_vals.end();
  int nread;
  char *line = NULL; size_t N=0;
  while (1) {
    if (getline(&line,&N,fid)==-1) break;
    nread = sscanf(line,"%lf %lf",&u[0],&u[1]);
    assert(nread==NDIM);
    node++;
    q = coup_vals.find(node);
    if (q!=qe) {
      int ext_node_indx = q->second;
      ext_node_data[ext_node_indx].u[0] = u[0];
      ext_node_data[ext_node_indx].u[1] = u[1];
    }
  }
  fclose(fid);

  // Computes potential on the external side, integrating
  // the tangential component of velocity
  int nnod_ext = ext_node_data.size();
  FastMat2 uu(1,NDIM), u1(1,NDIM), x1(1,NDIM), dx(1,NDIM), dpot_fm;
  ext_node_data[0].phi = 0.;
  for (int k=1; k<nnod_ext; k++) {
    u1.set(ext_node_data[k-1].x);
    uu.set(ext_node_data[k].x).add(u1).scale(0.5);

    x1.set(ext_node_data[k-1].u);
    dx.set(ext_node_data[k].x).rest(x1);
    
    double dpot = dpot_fm.prod(dx,uu,-1,-1).get();
    ext_node_data[k].phi = ext_node_data[k-1].phi + dpot;
  }

  for (int k=0; k<nnod_ext; k++) {
    ext_node &e = ext_node_data[k];
    printf("k=%d, vn %d, pn %d, x=%f %f, uu=%f %f, phi=%f\n",
	   k,e.vn,e.pn,e.x[0],e.x[1],e.u[0],e.u[1],e.phi);
  }
  PetscFinalize();
  exit(0);
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
  void init(TextHashTable *thash) { }
  double eval(double) { }
  ~coupling() { }
};

DEFINE_EXTENDED_AMPLITUDE_FUNCTION(coupling);
