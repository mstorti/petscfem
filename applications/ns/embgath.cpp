//__INSERT_LICENSE__
//$Id: embgath.cpp,v 1.26 2003/01/07 10:33:42 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/linkgraph.h>

#include "embgath.h"
extern Mesh *GLOBAL_MESH;
extern int MY_RANK,SIZE;

double my_power(double x,void* p) {
  return pow(x,*(double *)p);
}

class Cloud {
private:
  int nderiv, npol, nx;
  double nderiv_fact;
  FastMat2 A, xi, H, iH, AA;
public:
  void init(int nx_a, int nderiv_a,int npol_a) {
    nderiv = nderiv_a;
    npol = npol_a;
    nx = nx_a;
    A.resize(2,nx,npol+1);
    AA.resize(2,nx,npol+1);
    xi.resize(1,nx);
    H.resize(2,npol+1,npol+1);
    iH.resize(2,npol+1,npol+1);
    nderiv_fact = 1.;
    for (int k=2; k<=nderiv; k++)
      nderiv_fact *= double(k);
  }
  void coef(FastMat2 &x, FastMat2 &w,double x0=0.) {
    A.ir(2,1).set(1.);
    xi.set(x).add(-x0);
    double h = xi.max_abs_all();
    xi.scale(1./h);
    for (int k=1; k<=npol; k++) {
      double exp = double(k);
      A.ir(2,k+1).set(xi).fun(&my_power,&exp);
    }
    A.rs();
    H.prod(A,A,-1,1,-1,2);
    iH.inv(H);
    AA.prod(A,iH,1,-1,2,-1);
    AA.ir(2,nderiv+1);
    w.set(AA).scale(nderiv_fact/pow(h,nderiv));
    AA.rs();
  }
  void clear() { 
    A.clear(); 
    xi.clear(); 
    H.clear();
    iH.clear();
  }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Contains the three non-trivial face orientations. All
    other may be obtained by reflections ad rotations. */
const int Quad2Hexa::faces[][8] = {
  0,1,2,3,4,5,6,7,
  1,5,6,2,0,4,7,3,
  0,4,5,1,3,7,6,2};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Remaps the volume connectivity #vol_conn# so that the face
    is in a standard position. 
    @param surf_map (input) the nodes of the volume element
    that are on the face. 
    @param vol_conn (input/output) the connectivity of the
    volume element. Is remapped so that the face is in a
    standard position. 
*/ 
int Surf2Vol::map_mask(const int *surf_map,int *vol_conn) {
  int nel_surf, nel_vol, nf = nfaces(nel_surf,nel_vol);
  int match=0;
  const int *fc, *vol;
  // Loop over all face orientation posibilities until
  // one of them matches
  for (int f=0; f<nf; f++) {
    face(f,fc,vol);
    match = 1;
    for (int l=0; l<nel_surf; l++) {
      if (surf_map[l] != fc[l]) {
	match=0;
	break;
      }
    }
    if (match) break;
  }
  // Verify that one of the orientations must match
  if (!match) return 0;
  // Remap volume connectivity
  vector<int> vol_conn_c(nel_vol);
  for (int j=0; j<nel_vol; j++) vol_conn_c[j] = vol_conn[vol[j]];
  for (int j=0; j<nel_vol; j++) vol_conn[j] = vol_conn_c[j];
  return 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Quad2Hexa::face(int j,const int *&fc,const int *&vol_ret) { 
  // Changes quad numbering orientation. 
  int spin_map[] = {0,3,2,1};
  // Splits face number j into a triplet: spin(0/1) - rota (0-3) - k (0-2)
  int spin,rota,k,m;
  spin = modulo(j,2,&m);
  rota = modulo(m,4,&k);
  // Construct volume connectivity corresponding to k/rota
  for (int l=0; l<4; l++) {
    int ll = modulo(l+rota,4);
    vol[l] = faces[k][ll];
    vol[l+4] = faces[k][ll+4];
  }
  // Change orientation if needed. 
  if (spin) {
    for (int l=0; l<4; l++) {
      vol_r[l] = vol[4+spin_map[l]];
      vol_r[4+l] = vol[spin_map[l]];
    }
    vol_ret = vol_r;
  } else vol_ret = vol;

  // Face are the first 4 nodes.
  // Rotate depending on `use_exterior_normal'
  for (int k=0; k<4; k++) 
    this_face[k] = (use_exterior_normal() ? 
		    vol_ret[spin_map[k]] : vol_ret[k]);
  fc = this_face;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Line2Quad::face(int j,const int *&fc,const int *&vol) {
  static int fc_c[2], vol_c[4];
  static const int vol_cc[4] = {0, 1, 3, 2};
  for (int k=0; k<2; k++) fc_c[k] = (j+k) % 4;
  for (int k=0; k<4; k++) vol_c[k] = vol_cc[k]+j % 4;
  fc = fc_c;
  vol = vol_c;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Surf2Vol::Surf2Vol(const char *geom,int ndim,int nel,
		   int npg,int mat_version=GP_NEWMAT,
		   int use_exterior_normal_a=0) 
  : GPdata(geom,ndim,nel,npg,mat_version),
    use_exterior_normal_m(use_exterior_normal_a) {}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int embedded_gatherer::ask(const char *jobinfo,int &skip_elemset) {
  // Only accepts the `gather' jobinfo. 
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(gather);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void embedded_gatherer::initialize() {
  // For each face element, finds the volume element
  int ierr;
  //o The name of the volume elemset where to find
  // the volume element.
  TGETOPTDEF_S(thash,string,volume_elemset,);

  // Verifiy that the elemset name was given
  PETSCFEM_ASSERT0(volume_elemset.size()>0,
		   "embedded_gatherer: volume element name not given\n");
  int nelemsets = da_length(GLOBAL_MESH->elemsetlist);

  // Find volume elemset
  Elemset *vol_elem = GLOBAL_MESH->find(volume_elemset);
  // Verifiy that the elemset was found
  PETSCFEM_ASSERT(vol_elem,"Can't find volume element name: %s\n",
		  volume_elemset.c_str())

  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  // `npg_c' is a (dirty) trick to avoid collision between local
  // npg name and `npg' name in class :-(
  { int &npg_c = npg;
  //o Number of Gauss points.
  TGETOPTNDEF(thash,int,npg,none);
  npg_c = npg;
  }
  // ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  TGETOPTNDEF(thash,int,ndim,none); //nd
  //o Use exterior or interior normal
  TGETOPTDEF(thash,int,use_exterior_normal,1);
  //o Identify automatically the internal volume elements with a face
  // on the surface
  TGETOPTDEF(thash,int,identify_volume_elements,0);
  //o Number of layers in the normal direction.
  TGETOPTDEF_ND(thash,int,layers,1);
  PETSCFEM_ASSERT0(layers>=1,
		   "embedded_gatherer: Number of layers must be integer >=1\n");
  PETSCFEM_ASSERT(layers<=3,"embedded_gatherer: not supported yet layers>2,"
		  " entered layers: %d\n",layers);

  int ndimel=ndim-1;
  if (geometry=="quad2hexa") {
    sv_gp_data = new Quad2Hexa(geometry.c_str(),ndim,nel,npg,
			       GP_FASTMAT2,use_exterior_normal);
  } else if (geometry=="line2quad") {
    sv_gp_data = new Line2Quad(geometry.c_str(),ndim,nel,npg,
			       GP_FASTMAT2,use_exterior_normal);
  } else PETSCFEM_ERROR("embedded_gatherer: unknown geometry %s\n",geometry.c_str());

  surface_nodes(nel_surf,nel_vol);
  assert(nel_surf>0 && nel_surf<=nel);
  assert(nel_vol <= nel); //
  assert(nel_vol <= vol_elem->nel);
  assert(nel == nel_surf*(layers+1));
  if (identify_volume_elements) {
    assert(2*nel_surf==nel_vol);
    for (layer=0; layer<layers; layer++) {
      // Mark nodes on the surface
      int nnod = GLOBAL_MESH->nodedata->nnod;
      // surface:= is surface[k]==0 then k is not on the surface
      // if != 0 then surface[k] is the number of surface node +1
      vector<int> surface(nnod,0);
      // maps surface numbering (0 to surf_nodes-1) to global (0 to nnod-1)
      vector<int> srf2glb;
      for (int e=0; e<nelem; e++) {
	int *icorow = icone + nel*e + layer* nel_surf;
	for (int j=0; j<nel_surf; j++) surface[icorow[j]-1]=1;
      }
      // Count surface nodes
      int surf_nodes = 0;
      for (int k=0; k<nnod; k++) {
	if (surface[k]) {
	  surface[k] = ++surf_nodes;
	  srf2glb.push_back(surf_nodes);
	}
      }

      // Construct graph for volume elemset
      LinkGraph graph;
      graph.set_chunk_size(10000);
      graph.init(surf_nodes);

      // Construct node to element array for the volume elemset
      for (int e=0; e<vol_elem->nelem; e++) {
	int *icorow = vol_elem->icone + vol_elem->nel*e;
	for (int j=0; j<nel_vol; j++) {
	  int node = icorow[j]-1;
	  int snode = surface[node]-1;
	  if (snode>=0) graph.add(snode,e);
	}
      }

      // For each surface element look for the corresponding
      // volume element that shares a face
      vector<int> mask(nel_vol), icorow_c(nel_vol);
      for (int e=0; e<nelem; e++) {
	int *icorow = icone + nel*e + layer*nel_surf;
	LinkGraphRow row;
	assert(nel_surf>0);
	// Take list for first node
	int node = icorow[0]-1;
	int sf_node = surface[node]-1;
	graph.set_ngbrs(sf_node,row);
	LinkGraphRow::iterator q;
	int found=0;
	int *vicorow;
	int match=0;
	for (q=row.begin(); q!=row.end(); q++) {
	  int ve = *q; // the volume element
	  vicorow = vol_elem->icone + vol_elem->nel*ve;
	  for (int j=0; j<nel_vol; j++) mask[j]=-1;
	  found=0;
	  for (int j=0; j<nel_surf; j++) {
	    int sf_node = icorow[j];
	    for (int k=0; k<nel_vol; k++) {
	      if (vicorow[k]==sf_node) {
		mask[j] = k;
		break;
	      }
	    }
	    if (mask[j]==-1) break;
	    found++;
	  }
	  if (found==nel_surf) {
	    for (int j=0; j<nel_vol; j++) icorow_c[j] = vicorow[j];
	    // Volume element was found, find map and
	    if (sv_gp_data->map_mask(mask.begin(),icorow_c.begin())) {
	      match=1;
	      for (int j=0; j<nel_vol; j++) icorow[j] = icorow_c[j];
	      break;
	    }
	  }
	}
	if (!match) {
	  PetscPrintf(PETSC_COMM_WORLD,
		      "embedded_gatherer: Can't find matching volume element"
		      " to surface element %d\n",e);
	  PetscFinalize();
	  exit(0);
	}
      }

#if 0
      if (MY_RANK==0) {
	printf("Surface element connectivities: \n");
	for (int e=0; e<nelem; e++) {
	  int *icorow = icone + nel*e;
	  printf("surf.el. %d: ",e+1);
	  for (int j=0; j<nel_vol; j++) printf("%d ",icorow[j]);
	  printf("\n");
	}
      }
#endif
 
      graph.clear();
      surface.clear();
      srf2glb.clear();
    }
  }
}
  
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int embedded_gatherer::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
		      Dofmap *dofmap,const char *jobinfo,int myrank,
		      int el_start,int el_last,int iter_mode,
		      const TimeData *time) {

  int ierr;

  GET_JOBINFO_FLAG(gather);
  assert(gather);

  //o Position in gather vector
  TGETOPTDEF(thash,int,gather_pos,0);
  //o How many gather values will be contributed by this elemset
  TGETOPTDEF_ND(thash,int,gather_length,0);
  //o Dimension of the embedding space
  TGETOPTNDEF(thash,int,ndim,none);
  int ndimel = ndim-1;

#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 

#if 0
#define DSHAPEXI (*(*sv_gp_data).FM2_dshapexi[ipg])
#define SHAPE    (*(*sv_gp_data).FM2_shape[ipg])
#define WPG      ((*sv_gp_data).wpg[ipg])
#else
#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])
#endif

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)

  // get number of fields per node (constant+variables)
  int nu=nodedata->nu;

  int ja = 0;
  double *locst = arg_data_v[ja++].locst;
  double *locst2 = arg_data_v[ja++].locst;
  int options = arg_data_v[ja].options;
  vector<double> *values = arg_data_v[ja++].vector_assoc;
  int nvalues = values->size();
  assert(gather_pos+gather_length <= nvalues); // check that we don't put values
				   // beyond the end of global vector
				   // `values'
  vector<double> pg_values(gather_length);

  FastMat2 xloc(2,nel,ndim);

  Cloud cloud;
  cloud.init(layers+1,1,layers);
  
  FastMat2 Jaco(2,ndim,ndim),
    iJaco(2,ndim,ndim),staten(3,layers+1,nel_surf,ndof), 
    stateo(3,layers+1,nel_surf,ndof),
    u_old_l(2,layers+1,ndof),u_l(2,layers+1,ndof),
    u(1,ndof), u_old(1,ndof),
    n(1,ndim),xpgl(2,layers+1,ndim),xpg(1,ndim),grad_u(2,ndim,ndof),
    grad_uold(2,ndim,ndof),dshapex(2,ndim,nel),
    xn(1,layers+1),w(1,layers+1),state_pg(2,layers+1,ndof),
    grad_u_xi(2,ndim,ndof),grad_uold_xi(2,ndim,ndof);

  Time * time_c = (Time *)time;
  double t = time_c->time();
  
  // Initialize the call back functions
  init();

  GPdata gp_data("cartesian2d",ndim,4,npg,GP_FASTMAT2);
  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  int ielh=-1,kloc;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();
    ielh++;

    xloc.reshape(2,nel,ndim);
    for (kloc=0; kloc<nel; kloc++) {
      int node = ICONE(k,kloc);
      xloc.ir(1,kloc+1).set(&NODEDATA(node-1,0));
    }
    xloc.rs();
    xloc.reshape(3,layers+1,nel_surf,ndim);

    staten.set(&(LOCST(ielh,0,0)));
    stateo.set(&(LOCST2(ielh,0,0)));

    // Let user do some things when starting with an element
    element_hook(k);

    for (int ipg=0; ipg<npg; ipg++) {
      FastMat2 &shape = SHAPE;
      FastMat2 &dshapexi = DSHAPEXI;
      // Gauss point coordinates in several layers (layers+1 x npg x ndim)
      xpgl.prod(shape,xloc,-1,1,-1,2);

      // Jacobian master coordinates -> real coordinates (on surface)
      Jaco.is(1,1,2);
      xloc.ir(1,1);
      Jaco.prod(dshapexi,xloc,1,-1,-1,2);
      xloc.rs();

      // In surface jacobian
      double detJaco;
      detJaco = Jaco.detsur(&n);
      Jaco.rs();
      n.scale(1./detJaco);
      if (detJaco <= 0.) {
	printf("Jacobian of element %d is negative or null\n"
	       " Jacobian: %f\n",k,detJaco);
	PetscFinalize();
	exit(0);
      }
      double wpgdet = detJaco*WPG;

      // `w' is the coefficients to compute the derivative of a function with coordinates `xn' 
      xn.prod(n,xpgl,-1,1,-1);
      xn.add(-xn.get(1));
      cloud.coef(xn,w);
      
      // 3D Jacobian
      Jaco.ir(1,3).prod(xpgl,w,-1,1,-1).rs();
      iJaco.inv(Jaco);
      
      // Values and Gradients of variables at Gauss point
      // new state (t^n+1)
      staten.ir(1,1);
      u.prod(shape,staten,-1,-1,1);
      grad_u_xi.is(1,1,ndimel).prod(dshapexi,staten,1,-1,-1,2).rs();
      staten.rs();
      // old state (t^n)
      stateo.ir(1,1);
      u_old.prod(shape,stateo,-1,-1,1);
      grad_uold_xi.is(1,1,ndimel).prod(dshapexi,stateo,1,-1,-1,2).rs();
      stateo.rs();

      // state in all layers
      state_pg.prod(shape,staten,-1,1,-1,2);
      grad_u_xi.ir(1,ndim).prod(state_pg,w,-1,1,-1).rs();
      grad_u.prod(iJaco,grad_u_xi,1,-1,-1,2);
      
      state_pg.prod(shape,stateo,-1,1,-1,2);
      grad_uold_xi.ir(1,ndim).prod(state_pg,w,-1,1,-1).rs();
      grad_uold.prod(iJaco,grad_uold_xi,1,-1,-1,2);
      
      xpgl.ir(1,1); xpg.set(xpgl); xpgl.rs();
      set_pg_values(pg_values,u,u_old,grad_u,grad_uold,xpg,n,wpgdet,t);
      if (options & VECTOR_ADD) {
	for (int j=0; j<gather_length; j++) {
	  (*values)[gather_pos+j] += pg_values[j];
	}
      } else assert(0);
    }
  }  
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void visc_force_integrator::init() {
  int ierr;
  FastMat2 elc;
  //o Dimension of the embedding space
  TGETOPTNDEF(thash,int,ndim,none);
  ndim_m = ndim;

  compute_moment = (gather_length==2*ndim);
  force.resize(1,ndim);
  moment.resize(1,ndim);
  x_center.resize(1,ndim).set(0.);
  dx.resize(1,ndim);
  strain_rate.resize(2,ndim,ndim);
  sigma.resize(2,ndim,ndim);
  //o Viscosity of the fluid
  TGETOPTDEF_ND(thash,double,viscosity,0.);
  PETSCFEM_ASSERT0(viscosity>0.,"Viscosity should be positive.");  
  //o _T: double[ndim] _N: moment_center _D: null vector 
  // _DOC: Center of moments. _END
  get_double(thash,"moment_center",x_center.storage_begin(),1,ndim);  
  //o Rest the component of this rigid movement. 
  TGETOPTDEF_ND(thash,double,viscosity,0.);
  // Rotation angular velocity 
  Omega.resize(1,ndim).set(0.);
  ierr = get_double(thash,"Omega",Omega.storage_begin(),1,ndim);
  // Velocity gradient corresponding to the rigid movement
  elc.eps_LC();
  rigid_grad_u.prod(elc,Omega,1,2,-1,-1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void visc_force_integrator
::set_pg_values(vector<double> &pg_values,FastMat2 &u,
		FastMat2 &uold,FastMat2 &grad_u, FastMat2 &grad_uold, 
		FastMat2 &xpg,FastMat2 &n,
		double wpgdet,double time) {

  //#define SHV(pp) pp.print(#pp ": ")
#define SHV(pp) {}
  SHV(grad_u);
  grad_u.is(2,1,ndim_m).rest(rigid_grad_u);
  strain_rate.set(grad_u);
  grad_u.t();
  strain_rate.add(grad_u).scale(0.5);
  grad_u.rs();
  SHV(strain_rate);
  sigma.set(strain_rate).scale(2.*viscosity);
  sigma.d(1,2).add(-u.get(ndim_m+1)).rs();
  SHV(sigma);
  
  // Force contribution = normal * pressure * weight of Gauss point
  force.prod(sigma,n,1,-1,-1).scale(wpgdet);
  // export forces to return vector
  force.export_vals(pg_values.begin());
  SHV(force);

  if (compute_moment) {
    // Position offset of local point to center of moments
    SHV(xpg);
    SHV(x_center);
    dx.set(xpg).rest(x_center);
    SHV(dx);
    // Moment contribution = force X dx
    moment.cross(force,dx);
    SHV(moment);
    // export forces to return vector
    moment.export_vals(pg_values.begin()+ndim_m);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void visc_force_integrator
::surface_nodes(int &nel_surf,int &nel_vol) { 
  nel_surf=4; 
  nel_vol=8; 
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
