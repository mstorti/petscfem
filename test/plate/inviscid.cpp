//__INSERT_LICENSE__
//$Id: inviscid.cpp,v 1.8 2003/01/02 01:21:58 mstorti Exp $
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
  // viscous node, external potential node, external
  // potential node on second layer
  int vn, pn, pn1;
  // position, velocity, potential at node
  double x[NDIM], x2[NDIM], u[NDIM], phi, phi1;
};

double Uinf;
int inv_inc_lay;
typedef map<int,int> ext_map;
ext_map inv_indx_map,inv2_indx_map,visc_indx_map;
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
void coupling_visc_hook::time_step_pre(double t,int step) { 
  if (step>0) {
    printf("VISCOUS: waiting computed_step flag, step %d...\n",step);
    double step_sent = read_doubles(inv2visc,"computed_step");
    assert(int(step_sent)==step);
    printf("VISCOUS: received computed_step flag OK, step %d\n",step);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void coupling_visc_hook::time_step_post(double time,int step,
			       const vector<double> &gather_values) {
  fprintf(visc2inv,"computed_step %d\n",step,time);
}

DL_GENERIC_HOOK(coupling_visc_hook);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void coupling_inv_hook::init(Mesh &mesh,Dofmap &dofmap,
		     TextHashTableFilter *options_f,const char *name) {

  assert(SIZE==1); // Not implemented in parallel yet

  int ierr;
  TextHashTable *options = (TextHashTable *)TextHashTable::find(name);

  if (!MY_RANK) {
    printf("INVISCID: Opening fifos for communicating with VISCOUS.\n");

    visc2inv = fopen("visc2inv.fifo","r");
    assert(visc2inv);

    inv2visc = fopen("inv2visc.fifo","w");
    assert(inv2visc);
    setvbuf(inv2visc,NULL,_IOLBF,0);

    printf("INV2VISC: Done.\n");
    
    FILE *fid = fopen("ext.coupling_nodes.tmp","r");
    assert(fid);
    int nread, ext_node_indx=0;
    while (1) {
      int vn,pn;
      nread = fscanf(fid,"%d %d",&vn,&pn);
      if (nread!=2) break;
      ext_node_data.push_back(ext_node());
      ext_node_data[ext_node_indx].vn = vn;
      ext_node_data[ext_node_indx].pn = pn;
      inv_indx_map[pn] = ext_node_indx;      
      visc_indx_map[vn] = ext_node_indx;      
      ext_node_indx++;
    }
    assert(nread==EOF);
    fclose(fid);

    assert(fid = fopen("ext.con.tmp","r"));
    map<int, map<int, int> > layer_map;
    ext_map::iterator q,q1,qe;
    while (1) {
#define NEL 4
      int row[NEL];
      nread = fscanf(fid,"%d %d %d %d",&row[0],&row[1],&row[2],&row[3]);
      if (nread!=4) break;
      qe = inv_indx_map.end();
      for (int k=0; k<NEL; k++) {
	q =inv_indx_map.find(row[k]);
	if (q!=qe) {
	  for (int l=0; l<NEL; l++) {
	    q1 =inv_indx_map.find(row[l]);
	    if (q1==qe) {
	      layer_map[row[k]][row[l]]++;
	    }
	  }
	}
      }
    }
    assert(nread==EOF);
    fclose(fid);

    map<int,map <int,int> >::iterator r;
    map <int,int>::iterator s;
#if 1
    for (r=layer_map.begin(); r!=layer_map.end(); r++) {
      map <int,int> &m = r->second;
      for (s=m.begin(); s!=m.end(); s++) {
	printf("(%d,%d) -> %d\n",r->first,s->first,s->second);
      }
    }
#endif

    int nnod_ext = inv_indx_map.size();
    for (int k=1; k<nnod_ext-1; k++) {
      int pn = ext_node_data[k].pn;
      r = layer_map.find(pn);
      assert(r!=layer_map.end());
      int set_flag=0;
      map<int,int> &m = r->second;
      for (s=m.begin(); s!=m.end(); s++) {
	if (s->second==2) {
	  assert(!set_flag);
	  set_flag=1;
	  inv2_indx_map[s->first] = k;	  
	  ext_node_data[k].pn1 = s->first;
	}
      }
    }
    // First and last nodes
    for (int j=0; j<2; j++) {
      int k = (j==0? 0 : nnod_ext-1);
      int pn = ext_node_data[k].pn;
      r = layer_map.find(pn);
      assert(r!=layer_map.end());
      map<int,int> &m = r->second;
      assert(m.size()==2);
      int set_flag=0;
      for (s=m.begin(); s!=m.end(); s++) {
	int node = s->first;
	assert(s->second==1);
	if (inv2_indx_map.find(node)!=inv2_indx_map.end()) {
	  assert(!set_flag);
	  set_flag=1;
	  ext_node_data[k].pn1 = s->first;
	}
      }
    }

    for (int k=0; k<nnod_ext; k++) {
      ext_node &e = ext_node_data[k];
      printf("k %d, vn %d, pns %d %d\n",k,e.vn,e.pn,e.pn1);
    }
    PetscFinalize();
    exit(0);

    assert(fid = fopen("ext.nod.tmp","r"));
    double x[NDIM];
    int node=0;
    ext_map::iterator t, te=inv_indx_map.end();
    while (1) {
      nread = fscanf(fid,"%lf %lf",&x[0],&x[1]);
      if(nread!=NDIM) break;
      node++;
      t = inv_indx_map.find(node);
      if (t!=te) {
	ext_node_indx = t->second;
	printf("node %d, ext_node_indx %d, x %f %f\n",
	       node,ext_node_indx,x[0],x[1]);
	ext_node_data[ext_node_indx].x[0] = x[0];
	ext_node_data[ext_node_indx].x[1] = x[1];
      }
    }
    assert(nread==EOF);
    fclose(fid);
  }
  TGETOPTDEF_ND(options,double,Uinf,0.);
  assert(Uinf!=0.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void coupling_inv_hook::time_step_pre(double t,int step) { 
  printf("INVISCID: waiting computed_step flag, step %d...\n",step);
  double step_sent = read_doubles(visc2inv,"computed_step");
  assert(int(step_sent)==step);
  printf("INVISCID: received computed_step flag OK, step %d\n",step);
  FILE *fid = fopen("cylin.state.tmp","r");
  double u[NDIM];
  int node=0;
  ext_map::iterator q, qe=inv_indx_map.end();
  int nread;
  char *line = NULL; size_t N=0;
  while (1) {
    if (getline(&line,&N,fid)==-1) break;
    nread = sscanf(line,"%lf %lf",&u[0],&u[1]);
    assert(nread==NDIM);
    node++;
    q = visc_indx_map.find(node);
    if (q!=qe) {
      int ext_node_indx = q->second;
      ext_node_data[ext_node_indx].u[0] = u[0];
      ext_node_data[ext_node_indx].u[1] = u[1];
    }
  }
  fclose(fid);
  free(line); line=NULL; N=0;

  // Computes potential on the external side, integrating
  // the tangential component of velocity
  int nnod_ext = ext_node_data.size();
  FastMat2 uu(1,NDIM), u1(1,NDIM), x1(1,NDIM), dx(1,NDIM), dpot_fm;
  ext_node_data[0].phi = 0.;
  for (int k=1; k<nnod_ext; k++) {
    u1.set(ext_node_data[k-1].u);
    uu.set(ext_node_data[k].u).add(u1).scale(0.5);
    uu.addel(-Uinf,1);

    x1.set(ext_node_data[k-1].x);
    dx.set(ext_node_data[k].x).rest(x1);
    
    double dpot = dpot_fm.prod(dx,uu,-1,-1).get();
    ext_node_data[k].phi = ext_node_data[k-1].phi + dpot;
  }

#if 0
  for (int k=0; k<nnod_ext; k++) {
    ext_node &e = ext_node_data[k];
    printf("k=%d, vn %d, pn %d, x=%f %f, uu=%f %f, phi=%f\n",
	   k,e.vn,e.pn,e.x[0],e.x[1],e.u[0],e.u[1],e.phi);
  }
  PetscFinalize();
  exit(0);
#endif

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void coupling_inv_hook::time_step_post(double time,int step,
			       const vector<double> &gather_values) {
#if 0
  FILE *fid = fopen("ext.state.tmp","r");
  int node=0;
  ext_map::iterator q, qe=inv_indx_map.end();
  int nread;
  char *line = NULL; size_t N=0;
  while (1) {
    if (getline(&line,&N,fid)==-1) break;
    double phi1;
    nread = sscanf(line,"%lf",&phi1);
    assert(nread==1);
    node++;
    q = visc_indx_map.find(nodex-);
    if (q!=qe) {
      int ext_node_indx = q->second;
      ext_node_data[ext_node_indx].u[0] = u[0];
      ext_node_data[ext_node_indx].u[1] = u[1];
    }
  }
  fclose(fid);
  free(line); line=NULL; N=0;
#endif
  fprintf(inv2visc,"computed_step %d\n",step);
}

DL_GENERIC_HOOK(coupling_inv_hook);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class coupling : public DLGenericTmpl {
public:
  coupling() { }
  void init(TextHashTable *thash) { }
  double eval(double);
  ~coupling() { }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double coupling::eval(double) { 
  assert(field()==1);
  int node_c = node();
  ext_map::iterator q = inv_indx_map.find(node_c);
  assert(q!=inv_indx_map.end());
  int ext_node_indx = q->second;
  double phi = ext_node_data[ext_node_indx].phi;
  printf("node: %d, phi: %f\n",node_c,phi);
  return phi;
}

DEFINE_EXTENDED_AMPLITUDE_FUNCTION2(coupling);
