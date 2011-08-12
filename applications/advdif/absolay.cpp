//__INSERT_LICENSE__
// $Id: absolay.cpp,v 1.2.20.1 2007/02/23 19:18:07 rodrigop Exp $
#include "./absolay.h"

AbsorbingLayer *absorbing_layer_elemset_p = NULL;

void AbsorbingLayer::initialize() {
  nstep_histo = 4;
  ndof = 3;
}

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

  NSGETOPTDEF(int,ndim,0); //nd
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
  double alpha = 1.0;
  // The trapezoidal rule integration parameter
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
    alpha = (glob_param->alpha);

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

  //o Contains the Habso matrix: Hp (ndof x ndof),
  //  Hm (ndof x ndof),Uref (ndof)
  NGETOPTDEF_S(string,habso_file,"<none>");
  PETSCFEM_ASSERT0(habso_file!="<none>","habso_mat is required");  

  //o Scales volume damping term
  NSGETOPTDEF(double,Kabso,1.0);
  PETSCFEM_ASSERT0(Kabso>=0.0,
                   "Kabso must be non-negatvive");  
  
  // The case `ndimel=0' is a special case for concentrated
  // absorbing boundary condition, NOT distributed. 
  NSGETOPTDEF(int,ndimel,-1); //nd
  if (ndimel<0) ndimel = ndim;

  NSGETOPTDEF(int,use_layer,0); //nd
  // It seems that when using penalizatin you can't 
  // use alpha<1.0
  if (!use_layer) {
    alpha = 1.0;
    ndimel = 0;
  }

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
  FastMat2 Un(1,ndof),Uo(1,ndof),Ualpha(1,ndof),
    Jaco(2,ndimel,ndim),iJaco(2,ndimel,ndimel),
    dU(1,ndof),dshapex(2,ndimel,nel), normal(1,ndim),tmp1,tmp2,
    shape(1,nel);

  if (!flag) {
    flag = 1;
    Habso.resize(2,ndof,ndof);
    Uref.resize(1,ndof);
    dvector<double> habso_data;
    habso_data.cat(habso_file.c_str());
    PETSCFEM_ASSERT(habso_data.size()==2*ndof*ndof+ndof,
                    "habso_data must be 2*ndof*ndof. ndof %d, "
                    "habso_data size %d",ndof,habso_data.size());

    if (use_layer>0) Habso.set(habso_data.buff());
    else {
      Habso.set(habso_data.buff()+ndof*ndof);
    }

    Uref.set(habso_data.buff()+2*ndof*ndof);
    // FIXME:= this is not done if the master has not elements
    if (!MY_RANK) {
      Habso.print("Habso: ");
      Uref.print("Uref:");
    }
  }
  
  
#if 0
  Habso.print("Habso: ");
  PetscFinalize();
  exit(0);
#endif
  
  nen = nel*ndof;

  //o Type of element geometry to define Gauss Point data
  NGETOPTDEF_S(string,geometry,cartesian2d);
  GPdata gp_data;
  int npg=1;
  if (use_layer>0) {
    NSGETOPTDEF_ND(int,npg,0); //nd
    assert(npg>0);
    gp_data.init(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);
  }

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE	 (*gp_data.FM2_shape[ipg])
#define WPG	 (gp_data.wpg[ipg])

  double detJaco, wpgdet;
  int elem, ipg,node, jdim, kloc,lloc,ldof;

  // Position of current element in elemset and in chunk
  int k_elem, k_chunk;

  FastMatCacheList cache_list;
  if (use_fastmat2_cache) FastMat2::activate_cache(&cache_list);

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
    lstate.set(0.).axpy(lstaten,alpha).axpy(lstateo,(1-alpha));
    veccontr.set(0.);
    matloc.set(0.0);

    for (ipg=0; ipg<npg; ipg++) {
      
      if (use_layer>0) {
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
        shape.set(SHAPE);
      } else {
        shape.set(1.0);
        wpgdet = 1.0;
      }

      Un.prod(shape,lstate,-1,-1,1);
      Uo.prod(shape,lstateo,-1,-1,1);
      Ualpha.set(0.).axpy(Uo,1-alpha).axpy(Un,alpha);
      tmp1.prod(Habso,shape,shape,2,4,1,3);
      matloc.axpy(tmp1,alpha*Kabso*wpgdet);
      dU.set(Ualpha).minus(Uref);
      tmp2.prod(shape,Habso,dU,1,2,-1,-1);
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AbsorbingLayer::time_step_pre(int step) { 
  if (!MY_RANK) printf("in abso_hook::time_step_pre()\n");
  if (step==0) {
    PETSCFEM_ASSERT0(uhist.size()==0,"uhist should be empty here");  
    PETSCFEM_ASSERT0(whist.size()==0,"whist should be empty here");  
    u.cat("./STEPS/fsabso2d.state-0.tmp.gz");
    assert(u.size()%ndof==0);
    nnod = u.size()/ndof;
    if (!MY_RANK) printf("abso_hook: detected nnod %d\n",nnod);
    uhist.a_resize(3,nstep_histo,nnod,ndof);
    uhist.set(0.0);
    whist.a_resize(3,nstep_histo,nnod,ndof);
    whist.set(0.0);
  }
  char file[100];
  sprintf(file,"./STEPS/fsabso2d.state-%d.tmp.gz",step);
  u.read(file);
  PetscFinalize();
  exit(0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AbsorbingLayer::time_step_post() {
  if (!MY_RANK) printf("in abso_hook::time_step_post()\n");
}
