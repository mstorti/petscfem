//__INSERT_LICENSE__
//$Id merge-with-petsc-233-55-g52bd457 Fri Oct 26 13:57:07 2007 -0300$

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/linkgraph.h>
#include <src/cloud.h>
#include <src/surf2vol.h>
#include <src/surf2vol2.h>

#include "./embgath.h"
#include "./nsi_tet.h"

extern Mesh *GLOBAL_MESH;

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
  //o Position in gather vector
  TGETOPTDEF_ND(thash,int,gather_pos,0);
  //o How many gather values will be contributed by this elemset
  TGETOPTDEF_ND(thash,int,gather_length,0);
  pass_values_as_gather = gather_length>0;

  //o Store the computed values in per-element properties table
  TGETOPTDEF_ND(thash,int,pass_values_as_props,0);
  //o Number of values computed by the gatherer. 
  TGETOPTDEF(thash,int,store_values_length,gather_length);
  //o Name of property where gather values are stored
  TGETOPTDEF_S(thash,string,store_in_property_name,none);

  //o This option is relevant only if #pass_values_as_props#
  // is active. If #compute_densities==0# then the value set
  // in the per-element property is the integral of the
  // integrand over the element. Conversely, if
  // #compute_densities==1# the value set is the density of
  // the mean value of the integrand, i.e. the integral
  // divided by the area of the element. For instance, if
  // the integrand is the heat flow through the surface,
  // then if #compute_densities==0# then the value set in
  // the per-element property is the total heat flow through
  // the element (which has units of energy per unit time),
  // whereas if #compute_densities==1# then the value set is
  // the mean heat flow density (which has units of energy
  // per unit time and unit surface ). For the traction on a
  // surface, the passed value is the total force on the
  // element in one case (units of force) and the mean skin
  // friction in the other (units of force per unit
  // area). This option has no effect on the values passed
  // via the global #gather_values# vector. 
  TGETOPTDEF_ND(thash,int,compute_densities,0);

  nvalues = gather_length;
  if (pass_values_as_props) {
    if (!nvalues) nvalues = store_values_length;
    PETSCFEM_ASSERT0(store_in_property_name!="none",
                     "If `pass_values_as_props' is set, then "
                     "`store_in_property_name' is required!!");  
    props_hash_entry *phep;
    phep = (props_hash_entry *)
      g_hash_table_lookup(elem_prop_names,
                          store_in_property_name.c_str());
    PETSCFEM_ASSERT(phep," `store_in_property_name' was set "
                    "as %s, but no such property was defined in the "
                    "per element property table!",
                      store_in_property_name.c_str());
    phe = *phep;

    if (!nvalues) nvalues = phe.width;
    else PETSCFEM_ASSERT(phe.width==nvalues,
                         "length of property values does not match "
                         "the length required by gatherer\n"
                         "length[%s] = %d and nvalues= %d",
                         store_in_property_name.c_str(),phe.width,
                         nvalues);
    // MPI_Type_vector(nelem,nvalues,nelprops,MPI_DOUBLE,&stride);
    // MPI_Type_commit(&stride);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void embedded_gatherer
::before_assemble(arg_data_list &arg_datav,Nodedata *nodedata,
                  Dofmap *dofmap, const char *jobinfo,int myrank,
                  int el_start,int el_last,int iter_mode,
                  const TimeData *time_data) {
  if (pass_values_as_props) {
    // Clear props data corresponding to gather values
    for (int k=0; k<nelem; k++) {
      int l = k*nelprops+phe.position;
      for (int j=0; j<phe.width; j++) 
        elemprops[k+j] = 0.0;
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int embedded_gatherer::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
		      Dofmap *dofmap,const char *jobinfo,int myrank,
		      int el_start,int el_last,int iter_mode,
		      const TimeData *time) {

  int ierr;

  PETSCFEM_ASSERT0(!strcmp(jobinfo,"gather"),
                   "This routine is supposed to be called "
                   "only with `gather' jobinfo\n");  

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
  vector<double> *values = NULL;

  if (pass_values_as_gather) {
    values = arg_data_v[ja++].vector_assoc;
    
    // check that we don't put values beyond the end of global
    // vector `values'
    PETSCFEM_ASSERT(static_cast<unsigned int>(gather_pos+gather_length)
                    <=values->size(),
                    "Not enough positions in gather vector for this gatherer.\n"
                    "Element %s, gather_pos %d, gather_length %d\n"
                    "size of values vector (ngather) %s\n",
                    name(),gather_pos,gather_length,values->size());
  }

  vector<double> pg_values(nvalues);

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
  GPdata *gp_data_p=NULL;
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

    double area=0.0;
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
      if (detJaco<=0.) {
	detj_error(detJaco,k);
	set_error(1);
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
      area += wpgdet;
      set_pg_values(pg_values,u,u_old,grad_u,grad_uold,xpg,n,wpgdet,t);
      if (pass_values_as_gather) {
        if (options & VECTOR_ADD) {
          for (int j=0; j<nvalues; j++) {
            (*values)[gather_pos+j] += pg_values[j];
          }
        } else PETSCFEM_ERROR0("Doesn't make sense gather values "
                               "and !VECTOR_ADD");
      }
    }
    if (pass_values_as_props) {
      int l = k*nelprops+phe.position;
      for (int j=0; j<phe.width; j++) {
        elemprops[l+j] += pg_values[j];
        if (compute_densities) elemprops[l+j] /= area;
      }
    }
  }  
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  delete gp_data_p;
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void embedded_gatherer
::after_assemble(const char *jobinfo) {
  if (pass_values_as_props) {
    dvector<double> send,recv;
    send.mono(nvalues*nelem);
    send.reshape(2,nelem,nvalues);
    recv.mono(nvalues*nelem);
    recv.reshape(2,nelem,nvalues);
    for (int k=0; k<nelem; k++) {
      int l = k*nelprops+phe.position;
      for (int j=0; j<phe.width; j++)
        send.e(k,j) = elemprops[l+j];
    }
    //#define DBG
#ifdef DBG
    if (!MY_RANK) {
      printf("BEFORE ALLREDUCE\n");
      for (int k=0; k<nelem; k++) {
        printf("elem %d, vals ",k);
        int l = k*nelprops+phe.position;
        for (int j=0; j<phe.width; j++)
          printf(" %f",elemprops[l+j]);
        printf("\n");
      }
    }
#endif
    MPI_Allreduce(send.buff(),recv.buff(),send.size(),MPI_DOUBLE,
                  MPI_SUM,PETSCFEM_COMM_WORLD);
    for (int k=0; k<nelem; k++) {
      int l = k*nelprops+phe.position;
      for (int j=0; j<phe.width; j++)
        elemprops[l+j] = recv.e(k,j);
    }
    send.clear();
    recv.clear();
#ifdef DBG
    if (!MY_RANK) {
      printf("AFTER ALLREDUCE\n");
      for (int k=0; k<nelem; k++) {
        printf("elem %d, vals ",k);
        int l = k*nelprops+phe.position;
        for (int j=0; j<phe.width; j++)
          printf(" %f",elemprops[l+j]);
        printf("\n");
      }
    }
#endif
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void visc_force_integrator::init() {
  int ierr;
  FastMat2 elc;
  //o Dimension of the embedding space
  TGETOPTNDEF(thash,int,ndim,none);
  ndim_m = ndim;

  PETSCFEM_ASSERT0(ndim==2 || ndim==3,"Only for 2D or 3D");  
  if (ndim==3) {
    PETSCFEM_ASSERT(nvalues==3 || nvalues==6,
                    "In 3D nvalues should have length ndim "
                    "(compute forces only)\n"
                    "or 2*ndim (compute forces and moments)\n"
                    "nvalues = %d\n", nvalues);  
    compute_moment = (nvalues==6);
  } if (ndim==2) {
    PETSCFEM_ASSERT(nvalues==2 || nvalues==3,
                    "In 2D nvalues should have length 2 "
                    "(compute forces only)\n"
                    "or 3 (compute forces and moments)\n"
                    "nvalues = %d\n", nvalues);  
    compute_moment = (nvalues==3);
  }
  force.resize(1,ndim);
  if (ndim==3) moment.resize(1,ndim);
  else moment.resize(1,1);
  x_center.resize(1,ndim).set(0.);
  dx.resize(1,ndim);
  strain_rate.resize(2,ndim,ndim);
  sigma.resize(2,ndim,ndim);
  sigma_old.resize(2,ndim,ndim);
  //o Viscosity of the fluid
  TGETOPTDEF_ND(thash,double,viscosity,0.);
  //o Mask for deviatoric component of stress tensor 
  TGETOPTDEF_ND(thash,double,dev_comp_mask,1.);
  //o Mask for pressure (normal) component of stress tensor 
  TGETOPTDEF_ND(thash,double,pressure_comp_mask,1.);
  PETSCFEM_ASSERT0(viscosity>0.,"Viscosity should be positive.");  

  //o alpha for time rule integration
  TGETOPTDEF_ND(thash,double,alpha,1.0);

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

#undef SHV
#define SHV(pp) {}
  grad_uold.is(2,1,ndim_m).rest(rigid_grad_u);
  strain_rate.set(grad_uold);
  grad_uold.t();
  strain_rate.add(grad_uold).scale(0.5);
  grad_uold.rs();
  SHV(strain_rate);
  sigma_old.set(strain_rate).scale(2.*viscosity*dev_comp_mask);
  sigma_old.d(1,2).add(-uold.get(ndim_m+1)*pressure_comp_mask).rs();

  SHV(grad_u);
  grad_u.is(2,1,ndim_m).rest(rigid_grad_u);
  strain_rate.set(grad_u);
  grad_u.t();
  strain_rate.add(grad_u).scale(0.5);
  grad_u.rs();
  SHV(strain_rate);
  sigma.set(strain_rate).scale(2.*viscosity*dev_comp_mask);
  sigma.d(1,2).add(-u.get(ndim_m+1)*pressure_comp_mask).rs();
  SHV(sigma);
  
  // sigma at time t(n+alpha)
  sigma.scale(alpha).axpy(sigma_old,1.0-alpha).rs();

  //  printf("alpha at integrator = %f \n",alpha);

  // Force contribution = normal * pressure * weight of Gauss point
  force.prod(sigma,n,1,-1,-1).scale(wpgdet);
  // export forces to return vector
  force.export_vals(&*pg_values.begin());
  SHV(force);

  if (compute_moment) {
    // Position offset of local point to center of moments
    SHV(xpg);
    SHV(x_center);
    dx.set(xpg).rest(x_center);
    SHV(dx);
    // Moment contribution = force X dx
    moment.cross(dx,force);
    SHV(moment);
    // export forces to return vector
    moment.export_vals(&*pg_values.begin()+ndim_m);
  }
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
