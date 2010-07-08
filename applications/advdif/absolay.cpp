//__INSERT_LICENSE__
// $Id: absolay.cpp,v 1.2.20.1 2007/02/23 19:18:07 rodrigop Exp $
#include "./absolay.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "int advective::ask(char *,int &)"
int AbsorbingLayer::ask(const char *jobinfo,int &skip_elemset) {
   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_prof);
   return 0;
}

double rnd(double) {
  return 2.0*drand()-1.0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void AbsorbingLayer::
new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
	     const Dofmap *dofmap,const char *jobinfo,
	     const ElementList &elemlist,
	     const TimeData *time_data) try {

  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_prof);

  int ierr=0;

  int locdof,kldof,lldof;

  NSGETOPTDEF(int,npg,0); //nd
  NSGETOPTDEF(int,ndim,0); //nd
  assert(npg>0);
  assert(ndim>0);

  int nelprops,nel,ndof;
  elem_params(nel,ndof,nelprops);
  int nen = nel*ndof;

  // Unpack Dofmap
  int neq,nnod;
  neq = dofmap->neq;
  nnod = dofmap->nnod;

  // Unpack nodedata
  int nu=nodedata->nu;
  // H is a generalized local property passed per node with the nodal
  // coordinates. In shallow water nH =1 and H is the depth. It is
  // needed in order to compute the source term. In 1D Euler it may be
  // the area section of the tube. Its gradient is needed for the
  // source term in the momentum eqs.
  int nH = nu-ndim;
  FMatrix  Hloc(nel,nH),H(nH);

  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  // lambda_max:= the maximum eigenvalue of the jacobians.
  // used to compute the critical time step.
  vector<double> *dtmin;
  int jdtmin;
  GlobParam *glob_param=NULL;
  // The trapezoidal rule integration parameter
#define ALPHA (glob_param->alpha)
#define DT (glob_param->Dt)
  arg_data *staten=NULL,*stateo=NULL,*retval=NULL,
    *fdj_jac=NULL,*jac_prof=NULL,*Ajac=NULL;
  if (comp_res) {
    int j=-1;
    stateo = &arg_data_v[++j]; //[0]
    staten = &arg_data_v[++j]; //[1]
    retval  = &arg_data_v[++j];//[2]
    jdtmin = ++j;//[3]
#define DTMIN ((*(arg_data_v[jdtmin].vector_assoc))[0])
#define WAS_SET arg_data_v[jdtmin].was_set
    Ajac = &arg_data_v[++j];//[4]
    glob_param = (GlobParam *)arg_data_v[++j].user_data;;

    if (ADVDIF_CHECK_JAC)
      fdj_jac = &arg_data_v[++j];
  }

  FastMat2 matloc(4,nel,ndof,nel,ndof);

  //o Uses operations caches for computations with the FastMat2
  //  library for internal computations. Note that this affects also the
  //  use of caches in routines like fluxes, etc...
  NSGETOPTDEF(int,use_fastmat2_cache,1);

  //o Gravity
  NSGETOPTDEF(double,gravity,0.0);
  //o Reference water depth
  NSGETOPTDEF(double,h0,0.0);
  //o Scales volume damping term
  NSGETOPTDEF(double,Kabso,1.0);
  PETSCFEM_ASSERT0(Kabso>=0.0,
                   "Kabso must be non-negatvive");  

  // Initialize flux functions
  int ff_options=0;
  // adv_diff_ff->start_chunk(ff_options);
  // int ndimel = adv_diff_ff->dim();
  int ndimel = ndim;
  if (ndimel<0) ndimel = ndim;

  // adv_diff_ff->set_profile(prof_fields); // profile by equations
  matloc.set(1.0);
  if (comp_prof) {
    jac_prof = &arg_data_v[0];
    matloc.export_vals(jac_prof->profile);
  }

  int nlog_vars;
  const int *log_vars;
  // adv_diff_ff->get_log_vars(nlog_vars,log_vars);
  //o Use log-vars for $k$ and $\epsilon$
  NSGETOPTDEF(int,use_log_vars,0);
  if (!use_log_vars) nlog_vars=0;

  // Allocate local vecs
  FMatrix veccontr(nel,ndof),veccontr_mass(nel,ndof),
    xloc(nel,ndim),lstate(nel,ndof),
    lstateo(nel,ndof),lstaten(nel,ndof),dUloc_c(nel,ndof),
    dUloc(nel,ndof);
  FastMat2 Habso(2,ndof,ndof),Un(1,ndof),Uo(1,ndof),Ualpha(1,ndof),
    Uref(1,ndof),Jaco(2,ndimel,ndim),iJaco(2,ndimel,ndimel),
    dU(1,ndof),dshapex(2,ndimel,nel), normal(1,ndim),tmp1,tmp2;
#if 0
  static int flag=0;
  static FastMat2 Habsoc(2,ndof,ndof);
  if (!flag) {
    flag = 1;
    FastMat2 tmp3(2,ndof,ndof),tmp4(2,ndof,ndof);
    tmp3.fun(rnd);
    Habsoc.eye().axpy(tmp3,0.25);
    tmp4.ctr(tmp3,2,1);
    Habsoc.axpy(tmp4,0.25);
  }
  Habso.set(Habsoc);
#else
  Habso.eye();
#endif
  
#if 0
  Habso.print("Habso: ");
  PetscFinalize();
  exit(0);
#endif
  Uref.set(0.0).setel(h0,ndof);
  
  nen = nel*ndof;

  //o Type of element geometry to define Gauss Point data
  NGETOPTDEF_S(string,geometry,cartesian2d);
  GPdata gp_data(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE	 (*gp_data.FM2_shape[ipg])
#define WPG	 (gp_data.wpg[ipg])

  double detJaco, wpgdet;
  int elem, ipg,node, jdim, kloc,lloc,ldof;

  // Position of current element in elemset and in chunk
  int k_elem, k_chunk;

  FastMatCacheList cache_list;
  if (use_fastmat2_cache) FastMat2::activate_cache(&cache_list);

  // printf("[%d] %s start: %d last: %d\n",MY_RANK,jobinfo,el_start,el_last);
  for (ElementIterator element = elemlist.begin();
       element!=elemlist.end(); element++) try {

    element.position(k_elem,k_chunk);
    FastMat2::reset_cache();

    // Initialize element
    // adv_diff_ff->element_hook(element);

    // Get nodedata info (coords. etc...)
    element.node_data(nodedata,xloc.storage_begin(),
		      Hloc.storage_begin());

    if (comp_prof) {
      matloc.export_vals(element.ret_mat_values(*jac_prof));
      continue;
    }

    if (comp_res) {
      lstateo.set(element.vector_values(*stateo));
      lstaten.set(element.vector_values(*staten));
    }

    // State at time t_{n+\alpha}
    lstate.set(0.).axpy(lstaten,ALPHA).axpy(lstateo,(1-ALPHA));
    veccontr.set(0.);
    matloc.set(0.0);

    for (ipg=0; ipg<npg; ipg++) {
      
      //      Matrix xpg = SHAPE * xloc;
      Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);

      PETSCFEM_ASSERT0(ndim==ndimel,"not implemented yet");
      // iJaco.inv(Jaco);
      detJaco = Jaco.det();

      if (detJaco<=0.) {
        int k,ielh;
        element.position(k,ielh);
        detj_error(detJaco,k);
        set_error(1);
      }
      wpgdet = detJaco*WPG;
      // dshapex.prod(iJaco,DSHAPEXI,1,-1,-1,2);

      Un.prod(SHAPE,lstate,-1,-1,1);
      Uo.prod(SHAPE,lstateo,-1,-1,1);
      Ualpha.set(0.).axpy(Uo,1-ALPHA).axpy(Un,ALPHA);
      tmp1.prod(Habso,SHAPE,SHAPE,2,4,1,3);
      matloc.axpy(tmp1,ALPHA*Kabso*wpgdet);
      dU.set(Ualpha).rest(Uref);
      tmp2.prod(SHAPE,Habso,dU,1,2,-1,-1);
      veccontr.axpy(tmp2,-Kabso*wpgdet);
    }
    
    veccontr.export_vals(element.ret_vector_values(*retval));
    matloc.export_vals(element.ret_mat_values(*Ajac));
    
  } catch (GenericError e) {
    set_error(1);
    return;
  }

  FastMat2::void_cache();
  FastMat2::deactivate_cache();
} catch (GenericError e) {
  set_error(1);

}
