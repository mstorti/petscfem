//__INSERT_LICENSE__
//$Id: embgath.cpp,v 1.37 2003/02/27 03:32:41 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/linkgraph.h>
#include <src/cloud.h>
#include <src/surf2vol.h>
#include <src/surf2vol2.h>

#include "embgath.h"
extern Mesh *GLOBAL_MESH;
extern int MY_RANK,SIZE;

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
  } else if (geometry=="tri2prism") {
    sv_gp_data = new Tri2Prism(geometry.c_str(),ndim,nel,npg,
			       GP_FASTMAT2,use_exterior_normal);
  } else if (geometry=="line2quad") {
    sv_gp_data = new Line2Quad(geometry.c_str(),ndim,nel,npg,
			       GP_FASTMAT2,use_exterior_normal);
  } else PETSCFEM_ERROR("embedded_gatherer: unknown geometry \"%s\"\n",geometry.c_str());

  sv_gp_data->nfaces(nel_surf,nel_vol);
  assert(nel_surf>0 && nel_surf<=nel);
  assert(nel_vol <= nel); //
  assert(nel_vol <= vol_elem->nel);
  assert(nel == nel_surf*(layers+1));
  if (identify_volume_elements) 
    identify_volume_elements_fun(GLOBAL_MESH->nodedata->nnod,nel_surf,layers,nelem,
				 icone,nel,nel_vol,vol_elem,sv_gp_data);
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

  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,quad2hexa);
  GPdata *gp_data_p;
  if (geometry=="quad2hexa") 
    gp_data_p = new GPdata("cartesian2d",ndimel,4,npg,GP_FASTMAT2);
  else if (geometry=="tri2prism") 
    gp_data_p = new GPdata("triangle",ndimel,3,npg,GP_FASTMAT2);
  else if (geometry=="line2quad") 
    gp_data_p = new GPdata("cartesian1d",ndimel,2,npg,GP_FASTMAT2);
  else PETSCFEM_ERROR("Not known geometry %s for embgath",geometry.c_str());
  GPdata &gp_data = *gp_data_p;

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
      Jaco.is(1,1,ndimel);
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

      // `w' is the coefficients to compute the derivative of 
      //  a function with coordinates `xn' 
      xn.prod(n,xpgl,-1,1,-1);
      xn.add(-xn.get(1));
      cloud.coef(xn,w);
      
      // 3D Jacobian
      Jaco.ir(1,ndim).prod(xpgl,w,-1,1,-1).rs();
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
  delete gp_data_p;
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void visc_force_integrator::init() {
  int ierr;
  FastMat2 elc;
  //o Dimension of the embedding space
  TGETOPTNDEF(thash,int,ndim,none);
  ndim_m = ndim;

  assert(ndim==2 || ndim==3);
  if (ndim==3) {
    assert(gather_length==3 || gather_length==6);
    compute_moment = (gather_length==6);
  } if (ndim==2) {
    assert(gather_length==2 || gather_length==3);
    compute_moment = (gather_length==3);
  }
  force.resize(1,ndim);
  if (ndim==3) moment.resize(1,ndim);
  else moment.resize(1,1);
  x_center.resize(1,ndim).set(0.);
  dx.resize(1,ndim);
  strain_rate.resize(2,ndim,ndim);
  sigma.resize(2,ndim,ndim);
  //o Viscosity of the fluid
  TGETOPTDEF_ND(thash,double,viscosity,0.);
  PETSCFEM_ASSERT0(viscosity>0.,"Viscosity should be positive.");  
  if (0) {
    //o _T: double[ndim] _N: moment_center _D: null vector 
    // _DOC: Center of moments. _END
    get_double(thash,"moment_center",x_center.storage_begin(),1,ndim);  
    // Rotation angular velocity 
    Omega.resize(1,ndim).set(0.);
    ierr = get_double(thash,"Omega",Omega.storage_begin(),1,ndim);
    // Velocity gradient corresponding to the rigid movement
    elc.eps_LC();
    rigid_grad_u.prod(elc,Omega,1,2,-1,-1);
  } else {
    rigid_grad_u.resize(2,ndim,ndim).set(0.);
  }    
//    if (compute_moment) 
//      assert(pg_values.size() == (ndim==3? 6 ndim=2? 6 
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

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
