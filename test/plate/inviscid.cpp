//__INSERT_LICENSE__
//$Id: inviscid.cpp,v 1.12 2003/01/04 00:49:56 mstorti Exp $
#define _GNU_SOURCE

extern int MY_RANK,SIZE;

#include "./inviscid.h"
#include "./fifo.h"

extern int MY_RANK,SIZE;
#define NDIM 2    
#define NDOF (NDIM+1)    
//#define COUPLING_DEBUG

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class ext_node {
public:
  // viscous node, external potential node, external
  // potential node on second layer
  int vn, pn, pn1;
  // position, velocity, potential at node
  double x[NDIM], x1[NDIM], u[NDIM], phi, phi1;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class inv_ext_node {
public:
  int indx;
  double u_n;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double Uinf;
int inv_inc_lay;
typedef map<int,int> ext_map;
ext_map inv_indx_map,inv1_indx_map,visc_indx_map;
vector<ext_node> ext_node_data;
typedef map<int,inv_ext_node> coupling_visc_vel_t;
coupling_visc_vel_t coupling_visc_vel;
int computed_coupling_visc_vel=0;

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
  if (step>1) {
    PetscPrintf(PETSC_COMM_WORLD,
		"VISCOUS: waiting computed_step flag, step %d...\n",step);
    double step_sent = read_doubles(inv2visc,"computed_step");
    assert(int(step_sent)==step-1);
    PetscPrintf(PETSC_COMM_WORLD,
		"VISCOUS: received computed_step flag OK, step %d\n",step);

    FILE *fid = fopen("ext.coupling_normal_vel.tmp","r");
    assert(fid);
    FILE *fid2 = fopen("cylin.nod_fic_ext.tmp","r");
    assert(fid2);
    int nread;
    char *line = NULL; size_t N=0;
    int k=0;
    while (1) {
      if (getline(&line,&N,fid)==-1) break;
      int node; double u_n;
      nread = sscanf(line,"%d %lf",&node,&u_n);
      assert(nread==2);
      if (getline(&line,&N,fid2)==-1) break;
      nread = sscanf(line,"%d",&node);
      assert(nread==1);

      inv_ext_node & n = coupling_visc_vel[node];
      n.u_n = u_n;
      n.indx = k;
#ifdef COUPLING_DEBUG
      printf("VISCOUS: loads node %d, u_n %f, on coupling_visc_vel\n",
	     node,u_n);
#endif
      k++;
    }
    fclose(fid);
    fclose(fid2);
    free(line); line=NULL; N=0;
    computed_coupling_visc_vel=1;
  } else {
    PetscPrintf(PETSC_COMM_WORLD,
		"VISCOUS: step 1 do nothing ... Continue\n");
    computed_coupling_visc_vel=0;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void coupling_visc_hook::time_step_post(double time,int step,
			       const vector<double> &gather_values) {
  PetscPrintf(PETSC_COMM_WORLD,
	      "VISCOUS: sending computed_flag step %d\n",step);
  fprintf(visc2inv,"computed_step %d\n",step);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void coupling_visc_hook::close() {
  PetscPrintf(PETSC_COMM_WORLD,
	      "VISCOUS: sending computed_flag step -1 (stopping signal!!)\n");
  fprintf(visc2inv,"computed_step -1\n");
  fclose(visc2inv);
  read_doubles(inv2visc,"computed_step"); // read once floating value
  int step_sent = int(read_doubles(inv2visc,"computed_step"));
  printf("read from inv2visc: %d\n",step_sent);
  assert(step_sent==-1);
  fclose(inv2visc);
}

DL_GENERIC_HOOK(coupling_visc_hook);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void coupling_inv_hook::init(Mesh &mesh,Dofmap &dofmap,
			     TextHashTableFilter *options_f,
			     const char *name) {

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
#if 0
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
	  inv1_indx_map[s->first] = k;	  
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
	if (inv1_indx_map.find(node)!=inv1_indx_map.end()) {
	  assert(!set_flag);
	  set_flag=1;
	  ext_node_data[k].pn1 = s->first;
	}
      }
    }
    PetscPrintf(PETSC_COMM_WORLD,
		"INVISCID: identifies 2nd layer of nodes"
		" in inviscid region OK.\n");

#if 0
    for (int k=0; k<nnod_ext; k++) {
      ext_node &e = ext_node_data[k];
      printf("k %d, vn %d, pns %d %d\n",k,e.vn,e.pn,e.pn1);
    }
#endif

    assert(fid = fopen("ext.nod.tmp","r"));
    double x[NDIM];
    int node=0;
    ext_map::iterator t, te=inv_indx_map.end(),
      te1=inv1_indx_map.end();
    while (1) {
      nread = fscanf(fid,"%lf %lf",&x[0],&x[1]);
      if(nread!=NDIM) break;
      node++;
      t = inv_indx_map.find(node);
      if (t!=te) {
	ext_node_indx = t->second;
	ext_node_data[ext_node_indx].x[0] = x[0];
	ext_node_data[ext_node_indx].x[1] = x[1];
      }
      t = inv1_indx_map.find(node);
      if (t!=te1) {
	ext_node_indx = t->second;
	ext_node_data[ext_node_indx].x1[0] = x[0];
	ext_node_data[ext_node_indx].x1[1] = x[1];
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
  int step_sent = int(read_doubles(visc2inv,"computed_step"));
  if (step_sent==-1) {
    printf("INVISCID: received stopping signal... \n",step);
    fclose(visc2inv);
    printf("INVISCID: sending back acknowledge signal... \n",step);
    fprintf(inv2visc,"computed_step -1\n");
    fclose(inv2visc);
    PetscFinalize();
    exit(0);
  }
  assert(step_sent==step);
  printf("INVISCID: received computed_step flag OK, step %d\n",step);
  FILE *fid = fopen("cylin.state.tmp","r");
  assert(fid);
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
  FILE *fid = fopen("ext.state.tmp","r");
  assert(fid);
  int node=0;
  ext_map::iterator q, qe=inv_indx_map.end(),
    qe2=inv1_indx_map.end();
  int nread;
  char *line = NULL; size_t N=0;
  while (1) {
    if (getline(&line,&N,fid)==-1) break;
    double phi;
    nread = sscanf(line,"%lf",&phi);
    assert(nread==1);
    node++;
    q = inv_indx_map.find(node);
    if (q!=qe) {
      int ext_node_indx = q->second;
      ext_node_data[ext_node_indx].phi = phi;
    }
    q = inv1_indx_map.find(node);
    if (q!=qe2) {
      int ext_node_indx = q->second;
      ext_node_data[ext_node_indx].phi1 = phi;
    }
  }
  fclose(fid);
  free(line); line=NULL; N=0;

  int nnod_ext = ext_node_data.size();
  const int npoint=6, nterm=6;
  FastMat2 x(2,npoint,NDIM), xi(2,npoint,NDIM), psi(2,npoint,nterm), 
    jac(2,NDIM,NDIM), ijac(2,NDIM,NDIM), xc(1,NDIM), tmp(2,npoint,3),
    ipsi(2,npoint,nterm), a(1,nterm),
    pot(1,npoint),dpot_dxi(1,NDIM), dpot_dx(1,NDIM),
    t(1,NDIM), norm(1,NDIM), tmp1;
  double dphidn;

  assert(fid = fopen("ext.coupling_normal_vel.tmp","w"));
  for (int k=0; k<nnod_ext; k++) {
    int k_center = k + (k==0? 1 : k==nnod_ext-1? -1 : 0);
    xc.set(ext_node_data[k_center].x);

    x.ir(1,1).set(ext_node_data[k_center-1].x).rest(xc);
    x.ir(1,2).set(ext_node_data[k_center].x).rest(xc);
    x.ir(1,3).set(ext_node_data[k_center+1].x).rest(xc);
    x.ir(1,4).set(ext_node_data[k_center-1].x1).rest(xc);
    x.ir(1,5).set(ext_node_data[k_center].x1).rest(xc);
    x.ir(1,6).set(ext_node_data[k_center+1].x1).rest(xc);

    pot.setel(ext_node_data[k_center-1].phi,1);
    pot.setel(ext_node_data[k_center].phi,2);
    pot.setel(ext_node_data[k_center+1].phi,3);
    pot.setel(ext_node_data[k_center-1].phi1,4);
    pot.setel(ext_node_data[k_center].phi1,5);
    pot.setel(ext_node_data[k_center+1].phi1,6);

    x.ir(1,3);
    t.set(x);
    x.ir(1,1);
    t.rest(x);
    t.scale(1./t.norm_p_all(2.));
    norm.setel(-t.get(1),2);
    norm.setel(+t.get(2),1);

    x.ir(1,5);
    jac.ir(2,2).set(x);
    x.ir(1,3);
    jac.ir(2,1).set(x);
    x.ir(1,1);
    jac.rest(x).scale(0.5);
    x.rs();
    jac.rs();
    ijac.inv(jac);
    xi.prod(x,ijac,1,-1,2,-1);

    psi.ir(2,1).set(1.);
    xi.ir(2,1); psi.ir(2,2).set(xi);
    psi.ir(2,3).set(xi).mult(xi);
    psi.rs().is(2,1,3);
    tmp.set(psi);
    xi.ir(2,2);
    psi.rs().is(2,4,6).set(tmp);
    for (int l=4; l<=6; l++) psi.ir(2,l).mult(xi);
    psi.rs();
    xi.rs();

    ipsi.inv(psi);
    a.prod(ipsi,pot,1,-1,-1);
    dpot_dxi.setel(a.get(2),1).setel(a.get(4),2);
    dpot_dx.prod(ijac,dpot_dxi,-1,1,-1);
    // dpot_dx.print("dpot_dx: ");
    dphidn = tmp1.prod(dpot_dx,norm,-1,-1).get();

#ifdef COUPLING_DEBUG
    printf("node %d, indx %d, dphidn  %f\n",ext_node_data[k].vn,k,dphidn);
#endif
    fprintf(fid,"%d %f\n",ext_node_data[k].vn,dphidn);
  }
  fclose(fid);

  PetscPrintf(PETSC_COMM_WORLD,
	      "INVISCID: sending computed_flag step %d\n",step);
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
#ifdef COUPLING_DEBUG
  printf("node: %d, phi: %f\n",node_c,phi);
#endif
  return phi;
}

DEFINE_EXTENDED_AMPLITUDE_FUNCTION2(coupling);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class visc_coupling : public DLGenericTmpl {
public:
  visc_coupling() { }
  void init(TextHashTable *thash) { }
  double eval(double);
  ~visc_coupling() { }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double visc_coupling::eval(double) { 
  if (field()!=1 || !computed_coupling_visc_vel) return 0.;
  int node_c = node();
  coupling_visc_vel_t::iterator q = coupling_visc_vel.find(node_c);
  PETSCFEM_ASSERT(q != coupling_visc_vel.end(),
		  "Can't find node in coupling visc node list\n"
		  "node %d\n",node_c);
#ifdef COUPLING_DEBUG
  printf("VISCOUS: node %d, indx %d, normal vel. %f\n",
	 node_c,q->second.indx,q->second.u_n);
#endif
  return q->second.u_n;
}

DEFINE_EXTENDED_AMPLITUDE_FUNCTION2(visc_coupling);
