//__INSERT_LICENSE__
//$Id: genload.cpp,v 1.19 2002/08/27 02:53:46 mstorti Exp $
extern int comp_mat_each_time_step_g,
  consistent_supg_matrix_g,
  local_time_step_g;
extern int MY_RANK,SIZE;
 
#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/util2.h>
#include <src/fastmat2.h>

#include "advective.h"
#include "genload.h"

GenLoad::~GenLoad() {};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void HFilmFun::q(FastMat2 &uin,FastMat2 &flux,FastMat2 &jacin) {
  PETSCFEM_ERROR0("Not defined one layer flux function!\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void HFilmFun::q(FastMat2 &uin,FastMat2 &uout,FastMat2 &flux,
		 FastMat2 &jacin,FastMat2 &jacout) {
  PETSCFEM_ERROR0("Not defined two layers flux function!\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinearHFilmFun::element_hook(ElementIterator &element) {
  h->element_hook(element);
  s->element_hook(element);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int advective::ask(char *,int &)"
int GenLoad::ask(const char *jobinfo,int &skip_elemset) {
   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_prof);
   return 0;
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "double detsur(FastMat2 &Jaco, FastMat2 &S)"
double detsur(FastMat2 &Jaco, FastMat2 &S) {
  int n=Jaco.dim(2);
  int m=Jaco.dim(1);
  if (m==0) return 1.;
  static FastMat2 g;
  if (!g.is_defined()) g.resize(2,m,m);
  g.prod(Jaco,Jaco,1,-1,2,-1);
  return sqrt(g.det());
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
HFilmFun::HFilmFun(GenLoad *e) 
  : elemset(e), H_in(e->H_in), 
  H(e->H), H_out(e->H_out) {};


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "GenLoad::new_assemble"
void GenLoad::new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
			   const Dofmap *dofmap,const char *jobinfo,
			   const ElementList &elemlist,
			   const TimeData *time_data) {

  int ierr=0;

  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_prof);

  // nel2:= number of nodes in each side of the
  // element (if `double_layer' is in effect)
  int nelprops,nel,ndof,nel2;
  elem_params(nel,ndof,nelprops);
  int nen = nel*ndof;

  NSGETOPTDEF(int,npg,0); //nd
  NSGETOPTDEF(int,ndim,0); //nd
  //o The dimension of the element 
  NSGETOPTDEF(int,ndimel,ndim-1);
  PETSCFEM_ASSERT0(npg>0,"");
  PETSCFEM_ASSERT0(ndim>0,"");

  // nu:= total number of `constant' fields per node 
  // nH:= number of fields per node (not coordinates)
  // should be: nu = ndim + nH
  int nu=nodedata->nu;
  int nH = nu-ndim;
  if (nH>0) {
    H_m.resize(1,nH);
    H_out_m.resize(1,nH);
  }

  int locdof,kldof,lldof;
  char *value;

  //o Whether there is a double or single layer of nodes
  NSGETOPTDEF(int,double_layer,0);

  // allocate local vecs
  int kdof, kloc, node, jdim, ipg;
  FastMat2 veccontr(2,nel,ndof),xloc(2,nel,ndim),
    un(2,nel,ndof),uo(2,nel,ndof),ustar(2,nel,ndof),vecc2,
    h_in,h_out;
  FMatrix Hloc(nel,nH);

  nen = nel*ndof;
  FastMat2 matloc(4,nel,ndof,nel,ndof),matlocf(4,nel,ndof,nel,ndof),
    S(1,ndim),flux(1,ndof),load(1,ndof),jac_in,jac_out,Jaco,
    tmp1,tmp2,tmp3,tmp4;

  double detJ;
  FastMat2 u_in,u_out,U_in,U_out;

  int jdtmin;
  GlobParam *glob_param;
  // The trapezoidal rule integration parameter 
#define ALPHA (glob_param->alpha)
#define DT (glob_param->Dt)
  arg_data *staten,*stateo,*retval,*fdj_jac,*jac_prof,*Ajac;
  if (comp_res) {
    int j=-1;
    stateo = &arg_data_v[++j];
    staten = &arg_data_v[++j];
    retval  = &arg_data_v[++j];
    jdtmin = ++j;
#define DTMIN ((*(arg_data_v[jdtmin].vector_assoc))[0])
#define WAS_SET arg_data_v[jdtmin].was_set
    Ajac = &arg_data_v[++j];
    glob_param = (GlobParam *)arg_data_v[++j].user_data;;
#ifdef CHECK_JAC
    fdj_jac = &arg_data_v[++j];
#endif
  }

  if (comp_prof) {
    jac_prof = &arg_data_v[0];
    matlocf.set(1.);
  }

  h_film_fun->init(); // initialize hfilm function
  double hfilm;
  if (double_layer) {
    PETSCFEM_ASSERT0(nel % 2 ==0,"Number of nodes per element has to be even for "
		    "double_layer mode");
    nel2=nel/2;
    // state vector in both sides of the layer
    u_out.resize(2,nel2,ndof);
    jac_out.resize(2,ndof,ndof);
    U_out.resize(1,ndof);
  } else {
    nel2=nel;
  }
  // h_in:= h_out:= the `constant per node'
  // fields (not coordinates, neither unknowns) 
  if (nH>0) {
    h_in.resize(2,nel2,nH);
    h_out.resize(2,nel2,nH);
  }
  u_in.resize(2,nel2,ndof);
  jac_in.resize(2,ndof,ndof);
  U_in.resize(1,ndof);
  vecc2.resize(2,nel2,ndof);
  
  // Gauss Point data
  //o Type of element geometry to define Gauss Point data
  NGETOPTDEF_S(string,geometry,cartesian1d);
  GPdata gp_data(geometry.c_str(),ndimel,nel2,npg,GP_FASTMAT2);
  
#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);
  for (ElementIterator element = elemlist.begin(); 
       element!=elemlist.end(); element++) {

    FastMat2::reset_cache();

    // Load local node coordinates in local vector
    element.node_data(nodedata,xloc.storage_begin(),
		      Hloc.storage_begin());
    
    h_film_fun->element_hook(element); 

    if (comp_prof) {
      matlocf.export_vals(element.ret_mat_values(*jac_prof));
      continue;
    }

    matloc.set(0.);
    veccontr.set(0.);
    uo.set(element.vector_values(*stateo));
    un.set(element.vector_values(*staten));
    ustar.set(un).scale(ALPHA).axpy(uo,(1-ALPHA));
    
    ustar.is(1,1,nel2);
    u_in.set(ustar);
    ustar.rs();

    // Set h_in to the inner values of Hloc (if not double layer it is Hloc)
    if (nH>0) {
      Hloc.is(1,1,nel2);
      h_in.set(Hloc);
      Hloc.rs();
    }

    if (double_layer) {
      ustar.is(1,nel2+1,nel);
      u_out.set(ustar);
      ustar.rs();

      if (nH>0) {
	Hloc.is(1,nel2+1,nel);
	h_out.set(Hloc);
	Hloc.rs();
      }
    }

    // loop over Gauss points
    for (ipg=0; ipg<npg; ipg++) {

      xloc.is(1,1,nel2);
      Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);
      xloc.rs();

      detJ = Jaco.detsur();

      U_in.prod(SHAPE,u_in,-1,-1,1);
      if (nH>0) H_m.prod(SHAPE,h_in,-1,-1,1);
      if (double_layer) {
	U_out.prod(SHAPE,u_out,-1,-1,1);
	if (nH>0) H_out_m.prod(SHAPE,h_out,-1,-1,1);
	h_film_fun->q(U_in,U_out,flux,jac_in,jac_out);
      } else {
	h_film_fun->q(U_in,flux,jac_in);
      }
      
      double wpgdet = detJ*WPG;
      tmp1.set(SHAPE).scale(wpgdet);
      vecc2.prod(tmp1,flux,1,2);

      veccontr.is(1,1,nel2).add(vecc2);
      veccontr.rs();
      // Contribution to jacobian from interior side
      tmp2.set(SHAPE).scale(wpgdet);
      tmp3.prod(SHAPE,tmp2,1,2);
      tmp4.prod(tmp3,jac_in,1,3,2,4);
      matloc.is(1,1,nel2).is(3,1,nel2).add(tmp4);
      if (double_layer) {
	matloc.is(1).is(1,nel2+1,nel).rest(tmp4);
	veccontr.is(1,nel2+1,nel).rest(vecc2);
	veccontr.rs();
	// Contribution to jacobian from exterior side
	tmp4.prod(tmp3,jac_out,1,3,2,4);
	matloc.rs().is(1,1,nel2).is(3,nel2+1,nel).add(tmp4);
	matloc.is(1).is(1,nel2+1,nel).rest(tmp4);
      }
      matloc.rs();
    }

    veccontr.export_vals(element.ret_vector_values(*retval));
#ifdef CHECK_JAC
    veccontr.export_vals(element.ret_fdj_values(*fdj_jac));
#endif
    if (comp_mat_each_time_step_g) 
      matloc.export_vals(element.ret_mat_values(*Ajac));
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
