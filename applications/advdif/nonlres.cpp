/* $Id: nonlres.cpp,v 1.2.2.2 2004/02/25 14:25:42 mstorti Exp $ */

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/cloud.h>
#include <src/elemset.h>

#include "advective.h"
#include "nonlres.h"
#include "nwadvdif.h"

//extern Mesh *GLOBAL_MESH;

NonLinearRes::~NonLinearRes() {};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int NonLinearRes::ask(const char *jobinfo,int &skip_elemset) {
   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_prof);
   return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int LagrangeMult::ask(const char *jobinfo,int &skip_elemset) {
   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_prof);
   return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int NewAdvDif::ask(char *,int &)"
int AdvDiff_Abs_Nl_Res::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(comp_prof);
  DONT_SKIP_JOBINFO(comp_res);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void NewAdvDif::assemble"
void AdvDiff_Abs_Nl_Res::new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
				    const Dofmap *dofmap,const char *jobinfo,
				    const ElementList &elemlist,
				    const TimeData *time_data) {
  int nelprops,node,dof,ierr=0;
  elem_params(nel,ndof,nelprops);
  GET_JOBINFO_FLAG(comp_prof);
  GET_JOBINFO_FLAG(comp_res);
  NSGETOPTDEF(int,ndim,0);
  assert(ndim>0);
  int nu=nodedata->nu;
  int nH = nu-ndim;
  int ret_options=0;
  //  adv_diff_ff->start_chunk(ret_options); //ff ini

  //aqui en advdif esta primero en el arg_data_v el estado en t_n y despues en t_n+1.
  arg_data *locstold,*locst,*retval,*fdj_jac,*jac_prof,*Ajac;
  GlobParam *glob_param;
  double *hmin,Dt,rec_Dt;
  int ja_hmin;
#define WAS_SET arg_data_v[ja_hmin].was_set
  if (comp_res) {
    int ja=-1;
    locstold = &arg_data_v[++ja];
    locst = &arg_data_v[++ja];
    retval = &arg_data_v[++ja];
    ++ja; // Not used. Should be used for computation of automatic Dt
    Ajac = &arg_data_v[++ja];
    glob_param = (GlobParam *)(arg_data_v[++ja].user_data);
#ifdef CHECK_JAC
    fdj_jac = &arg_data_v[++ja];
#endif
    rec_Dt = 1./glob_param->Dt;
    if (glob_param->steady) rec_Dt=0.0;
  } else if (comp_prof) {
    jac_prof = &arg_data_v[0];
  }
  //o Using Lagrange multipliers leads to diagonal null terms, which can
  // cause zero pivots when using direct methods. With this option
  // a small term is added to the diagonal in order to fix this. The
  // term is added only in the Jacobian or also in the residual (which
  // results would be non-consistent). See option
  // \verb+lagrange_residual_factor+. 
  NSGETOPTDEF(double,lagrange_diagonal_factor,0.0);
  //o The diagonal term proportional to  \verb+lagrange_diagonal_factor+ 
  // may be also entered in the residual. If this is so
  // (\verb+lagrange_residual_factor=1+, then the
  // method is ``non-consistent'', i.e. the restriction is not exactly
  // satisfied by the non-linear scheme is exactly Newton-Raphson. If
  // not (\verb+lagrange_residual_factor=0+) then the restriction is
  // consistently satisfied but with a non exact Newton-Raphson. 
  NSGETOPTDEF(double,lagrange_residual_factor,0.0);
  //o Using Lagrange multipliers can lead to bad conditioning, which
  // causes poor convergence with iterative methods or amplification
  // of rounding errors. This factor scales the columns in the matrix
  // that correspond to the lagrange multipliers and can help in
  // better conditioning the system. 
  NSGETOPTDEF(double,lagrange_scale_factor,1.);
  FastMat2 matloc_prof(4,nel,ndof,nel,ndof),
    matloc(4,nel,ndof,nel,ndof), U(2,nel,ndof),R(2,nel,ndof);
  //o U is (U_{N},U_{N-1},..,U_{lagdof} U_{ref})
  if (comp_prof) matloc_prof.set(1.);
  nr = nres();
  //  AdvDiff_Abs_Nl_Res::init();
  FastMat2 r(1,nr),lambda(3,nel-2,ndof,nr),jac(3,nr,ndof,nel);
  jac.set(0.0);
  r.set(0.0);
  lambda.set(0.0);
  U.set(0.0); U_innodes.resize(2,nel-2,ndof).set(0.0);
  U_lagmul.resize(1,ndof).set(0.0);
  //////esto estaba en init()
  RI_.resize(2,nr,nel).set(0.);
  C_U_.resize(2,nr,nel).set(0.);
  drdU_.resize(3,nr,ndof,nel).set(0.);
  RI_ref.resize(1,ndof);
  C_U_ref.resize(1,ndof);
  drdU_ref.resize(2,nr,ndof);  
  U_ref.resize(1,ndof);
  xpe.resize(1,nel-3).set(0.);
  cpe.resize(1,nel-3).set(0.);
  RI_tmp.resize(1,nr).set(0.);
  drdU_tmp.resize(2,nr,ndof).set(0.);
  C_U_tmp.resize(1,nr).set(0.); //un rango menor que las primitivas
  int pp=adv_diff_ff->dim();
  normaln.resize(1,pp);
  adv_diff_ff->start_chunk(ret_options); //ff ini
  extr_cloud.init(nel-3,0,nel-4);
  double xne=0.;
  for (int i=1;i<=nel-3;i++) {
    xne=-1.*i;
    xpe.setel(xne,i); 
  }
  extr_cloud.coef(xpe,cpe);
  //////
#define COMPUTE_FD_RES_JACOBIAN
#ifdef COMPUTE_FD_RES_JACOBIAN
  FastMat2 res_fd_jac(3,nr,ndof,nel),res_pert(1,nr),U_pert(2,nel,ndof),
    lambda_pert(3,nel-2,ndof,nr),_fd_jac(3,nr,ndof,nel);
  res_fd_jac.set(0.0);res_pert.set(0.0);U_pert.set(0.0);lambda_pert.set(0.0);
  _fd_jac.set(0.0);
#endif
  
  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);
  for (ElementIterator element = elemlist.begin(); 
       element!=elemlist.end(); element++) {
    FastMat2::reset_cache();
    if(comp_prof) {
      matloc_prof.export_vals(jac_prof->profile);
      continue;
    }
    U.set(element.vector_values(*locst)); 
    matloc.set(0.0);
    R.set(0.0);
    if (comp_res) {
      init();
      element_hook(element);
      res(element,U,r,lambda,jac);
      U_innodes.set(U.is(1,1,nel-2));
      U.rs();
      U_lagmul.set(U.ir(1,nel-1));
      U.rs();
      R.is(1,1,nel-2).prod(lambda,U_lagmul,1,2,-1,-1).scale(lagrange_scale_factor);
      R.rs().ir(1,nel-1).is(2,1,nr).set(r)
	.axpy(U_lagmul,-lagrange_diagonal_factor*lagrange_residual_factor).rs();
      matloc.is(1,1,nel-2).ir(3,nel-1).is(4,1,nr).set(lambda)
	.scale(-lagrange_scale_factor).rs();
      jac.is(1,1,nr).is(2,1,ndof).is(3,1,nel-2);
      matloc.ir(1,nel-1).is(3,1,nel-2).is(2,1,nr).ctr(jac,1,3,2).scale(-1.).rs();
      matloc.ir(1,nel-1).ir(3,nel-1).d(2,4).is(2,1,nr)
	.set(lagrange_diagonal_factor).rs();
      jac.rs();
      R.export_vals(element.ret_vector_values(*retval));
      matloc.export_vals(element.ret_mat_values(*Ajac));
#ifdef CHECK_JAC
      R.export_vals(element.ret_fdj_values(*fdj_jac));
#endif
    }

#ifdef COMPUTE_FD_RES_JACOBIAN
    double eps_fd=1e-4;
    for (int jele=1; jele<=nel; jele++) {
      for (int jdof=1; jdof<=ndof; jdof++) {
	U_pert.set(U);      	
	U_pert.addel(eps_fd,jele,jdof);
	res(element,U_pert,res_pert,lambda_pert,_fd_jac);
	res_pert.rest(r).scale(1./eps_fd);
	res_fd_jac.ir(3,jele).ir(2,jdof).set(res_pert).rs();
      }
    }
#endif	
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  extr_cloud.clear();
}

void AdvDiff_Abs_Nl_Res::init() {
  get_prop(normaln_prop,"normaln");
  //  get_prop(U_ref_prop,"U_ref");
  assert(nel > 3);
  int ret_options=0;
  /*
    NSGETOPTDEF(string,vol_elemset,"none");
    assert(vol_elemset.length()>0);
    if (vol_elemset != "streamsw1d" || vol_elemset != "stream" || \
    vol_elemset != "advdif_swfm2t") {
    PetscPrintf(PETSC_COMM_WORLD,
    "Invalid value for \"volume_elemset\" option\n"
    "vol_elemset=\"%s\"\n",vol_elemset.c_str());
    PetscFinalize();
    exit(0);
    }
    Elemset *dummy_vol_elemset;
    dummy_vol_elemset = GLOBAL_MESH->find(vol_elemset);
    //check if found
    PETSCFEM_ASSERT(dummy_vol_elemset,"Can't find volume element name: %s\n",
    vol_elemset.c_str());
    // dynamic_cast from Elemset to streamsw1d (NewElemset)
    elemset_vol = dynamic_cast<const streamsw1d *>(dummy_vol_elemset);
    delete dummy_vol_elemset;
  */
}

void AdvDiff_Abs_Nl_Res::lag_mul_dof(int jr,int &node,int &dof) {
  //esto es por ahora
  node = 4;dof = jr;//creo que falta declararlas en la clase!!
}

void AdvDiff_Abs_Nl_Res::element_hook(ElementIterator &element){
  element_m = element;
  /* 
     assert(U_ref_prop.length==ndof);
     U_ref.set(prop_array(element_m,U_ref_prop));
  */
  int pp=adv_diff_ff->dim();
  assert(normaln_prop.length == pp);
  normaln.set(prop_array(element_m,normaln_prop));
}

void AdvDiff_Abs_Nl_Res::res(ElementIterator &element, FastMat2 &U, FastMat2 &r,
			   FastMat2 &lambda, FastMat2 &jac) {
  lambda.ir(1,1).eye();
  lambda.rs();
  for (int j=1;j<=nel;j++) {
    U.ir(1,j);
    adv_diff_ff->Riemann_Inv(U,normaln,RI_tmp,drdU_tmp,C_U_tmp);
    if (j == nel) {
      U_ref.set(U);
      RI_ref.set(RI_tmp);
      drdU_ref.set(drdU_tmp);
      C_U_ref.set(C_U_tmp);
    }
    U.rs();
    RI_.ir(2,j).set(RI_tmp).rs();
    C_U_.ir(2,j).set(C_U_tmp).rs();
    drdU_.ir(3,j).set(drdU_tmp).rs();
  }
  jac.set(drdU_);
  double c, rj=0., tmp=0., tmp2=0.;  
  C_U_.ir(2,1);
  for (int j=1;j<=nr;j++) {
    c = C_U_.get(j);
    assert(c != 0.0);
    FastMat2::branch();
    if (c > 0.0){
      FastMat2::choose(0);
      RI_.ir(2,j);
      rj = RI_.get(j);
      RI_.rs();
      RI_.ir(1,j);
      for(int k=1;k<nel-2;k++) {
	tmp = cpe.get(k); tmp2 = RI_.get(k+1);
	rj -= tmp*tmp2;
	drdU_.ir(1,j).ir(3,k+1);
       	jac.ir(1,j).ir(3,k+1).set(drdU_).scale(-tmp).rs();
	drdU_.rs();
      }
      jac.ir(1,j).ir(3,nel-1).set(0.).rs();
      RI_.rs();
      r.setel(rj,j);
    } else {
      FastMat2::choose(1);
      RI_.ir(2,1);
      tmp = RI_.get(j);tmp2 = RI_ref.get(j);
      rj = tmp-tmp2;
      RI_.rs();
      r.setel(rj,j);
      for(int k=1;k<=nel-1;k++) {
	jac.ir(1,j).ir(3,k+1);
	jac.set(0.).rs();
      }
    }
    FastMat2::leave();
    RI_.rs();
  }
  C_U_.rs();
  jac.rs();
}
