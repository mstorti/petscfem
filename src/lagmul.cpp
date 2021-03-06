//__INSERT_LICENSE__
/* $Id: lagmul.cpp,v 1.15.10.1 2007/02/19 20:23:56 mstorti Exp $ */

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <src/lagmul.h>

#define LagrangeMult GLagrangeMult

int ADVABSO_CURRENT_ELEMENT=INT_MAX;
extern TextHashTable *GLOBAL_OPTIONS;

LagrangeMult::~LagrangeMult() {}

void 
read_double_array_options(NewElemset &elemset,
			  const char *key,
			  vector<double> &v) {
  const char *line;
  elemset.get_entry(key,line);
  if(line) {
    v.clear();
    read_double_array(v,line);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// modif nsi_tet
#undef __FUNC__
#define __FUNC__ "LagrangeMult::new_assemble"
void LagrangeMult::
new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
	     const Dofmap *dofmap,const char *jobinfo,
	     const ElementList &elemlist,
	     const TimeData *time_data) try {

  nodedata_m = nodedata;
  int comp_mat, comp_mat_res;
  get_comp_flags(jobinfo,comp_mat,comp_mat_res);

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retvalmat,iele,j,nel,k,ndof,p,nel,q,ndof)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 

  // nr:= number of restrictions
  int nr,ierr=0,jr,jfic,dofic;
  // PetscPrintf(PETSCFEM_COMM_WORLD,"entrando a nsikeps\n");

  double lambda,rr;

  if (comp_mat_res) 
    get_data(arg_data_v,stateo,staten,retval,retvalmat);
  else if (comp_mat) 
    get_data(arg_data_v,retvalmat);

  time_m = double(* (const Time *) time_data);

  int nelprops,ndof;
  elem_params(nel,ndof,nelprops);

  vector<double> lagrange_diagonal_factor;
  //o _T: vector<double>
  // _N: lagrange_diagonal_factor
  // _D: scalar, 1e-5
  // _DOC: Using Lagrange multipliers leads to diagonal null
  // terms, which can cause zero pivots when using direct
  // methods. With this option a small term is added to the
  // diagonal in order to fix this. The term is added only
  // in the Jacobian or also in the residual (which results
  // would be non-consistent). See option
  // #lagrange_residual_factor# .
  //  _END
  read_double_array_options(*this,
			    "lagrange_diagonal_factor",
			    lagrange_diagonal_factor);

  //o The diagonal term proportional to   #lagrange_diagonal_factor#  
  // may be also entered in the residual. If this is so
  // ( #lagrange_residual_factor=1# , then the
  // method is ``non-consistent'', i.e. the restriction is not exactly
  // satisfied by the non-linear scheme is exactly Newton-Raphson. If
  // not ( #lagrange_residual_factor=0# ) then the restriction is
  // consistently satisfied but with a non exact Newton-Raphson. 
  NSGETOPTDEF(double,lagrange_residual_factor,0.);
  //o Using Lagrange multipliers can lead to bad conditioning, which
  // causes poor convergence with iterative methods or amplification
  // of rounding errors. This factor scales the columns in the matrix
  // that correspond to the lagrange multipliers and can help in
  // better conditioning the system. 
  NSGETOPTDEF(double,lagrange_scale_factor,1.);
  //o Using Lagrange multipliers can lead to bad conditioning, which
  // causes poor convergence with iterative methods or amplification
  // of rounding errors. This factor scales the row in the matrix
  // that correspond to the new equation. 
  NSGETOPTDEF(double,lagrange_row_scale_factor,1.);
  // Dimension of the embedding space (position vector of nodes)
  NSGETOPTDEF(int,ndim,0); //nd
  // Use or not caches for the FastMat2 libray
  NSGETOPTDEF(int,use_fastmat2_cache,1);

  // Call callback function defined by user initializing the elemset
  init();
  nr = nres(); // Get dimensions of problem

#if 0
  int nldf = lagrange_diagonal_factor.size();
  assert(nldf==0 || nldf==1 || nldf==nr);
  if (nldf==0) 
    lagrange_diagonal_factor.resize(nr,1e-5);
  if (nldf==1)
    lagrange_diagonal_factor
      .resize(nr,lagrange_diagonal_factor[0]);
  { static int flag=0;
  if (!flag) {
    flag=1;
    for (int j=0; j<nr; j++)
      printf("ldf[%d]=%g\n",j+1,lagrange_diagonal_factor[j]);
  }
  }
#endif
  
  FastMat2 matloc_prof(4,nel,ndof,nel,ndof),
    matloc(4,nel,ndof,nel,ndof), U(2,nel,ndof),R(2,nel,ndof),
    ldf_user(1,nr);
  xloc_m.resize(2,nel,ndim);
  if (comp_mat) matloc_prof.set(1.);

  // int nu=nodedata->nu;
  //int nH = nu-ndim;

  FastMat2 r(1,nr),w(3,nel,ndof,nr),jac(3,nr,nel,ndof);
  jac.set(0.);

  //#define COMPUTE_FD_RES_JACOBIAN
#ifdef COMPUTE_FD_RES_JACOBIAN
  FastMat2 res_fd_jac(3,nr,nel,ndof),
    d_res_fd_jac(3,nr,nel,ndof),
    res_pert(1,nr),U_pert(2,nel,ndof),
    lambda_pert(3,nel,ndof,nr),fd_jac(3,nr,nel,ndof);
  res_fd_jac.set(0.0);
  res_pert.set(0.0);
  U_pert.set(0.0);
  lambda_pert.set(0.0);
  fd_jac.set(0.0);
#endif
 
  FastMatCacheList cache_list;
  if (use_fastmat2_cache)
    FastMat2::activate_cache(&cache_list);
  int ielh=-1;
  nu=nodedata->nu;

  int k_chunk;
  for (element = elemlist.begin();
       element!=elemlist.end(); element++) try {

    element.position(elem,k_chunk);
    ADVABSO_CURRENT_ELEMENT = elem;
    FastMat2::reset_cache();
    ielh++;

    element_hook(element);
    if(comp_mat) {
      matloc_prof.export_vals(retvalmat->profile);
      continue;
    }      

    U.set(element.vector_values(*staten));
    matloc.set(0.);
    R.set(0.);

    if (comp_mat_res) {
      res(elem,U,r,w,jac);
      if (lagrange_diagonal_factor.size()==0)
	ldf_user.set(1e-5);
      else if (lagrange_diagonal_factor.size()==1)
	ldf_user.set(lagrange_diagonal_factor[0]);
      else if (lagrange_diagonal_factor.size()==(unsigned int)nr)
	ldf_user.set(&lagrange_diagonal_factor[0]);
      set_ldf(ldf_user,lagrange_diagonal_factor);
      for (jr=1; jr<=nr; jr++) {
	double ldf = ldf_user.get(jr);
	// get node/field of the Lag.mul.
	lag_mul_dof(jr,jfic,dofic);

	lambda = U.get(jfic,dofic);
	w.rs().ir(3,jr);
	R.axpy(w,lambda*lagrange_scale_factor);
	rr = r.get(jr) - ldf * 
	  lagrange_residual_factor*U.get(jfic,dofic);
	R.addel(rr*lagrange_row_scale_factor,jfic,dofic);

	matloc.rs().ir(3,jfic).
	  ir(4,dofic).axpy(w,-lagrange_scale_factor);
	matloc.rs().addel(ldf * lagrange_row_scale_factor,
			  jfic,dofic,jfic,dofic);
	jac.rs().ir(1,jr);
	matloc.ir(1,jfic).ir(2,dofic)
	  .axpy(jac,-lagrange_row_scale_factor);
      }
      R.rs().export_vals(element.ret_vector_values(*retval));
      matloc.rs().export_vals(element.ret_mat_values(*retvalmat));
      jac.rs();
      w.rs();
    }

#ifdef COMPUTE_FD_RES_JACOBIAN
    double eps_fd=1e-4;
    for (int jele=1; jele<=nel; jele++) {
      for (int jdof=1; jdof<=ndof; jdof++) {
	U_pert.set(U);      	
	U_pert.addel(eps_fd,jele,jdof);
	res(elem,U_pert,res_pert,lambda_pert,fd_jac);
	res_pert.minus(r).scale(1./eps_fd);
	res_fd_jac.ir(2,jele).ir(3,jdof)
	  .set(res_pert).rs();
      }
    }
    d_res_fd_jac
      .set(jac).minus(res_fd_jac);
    double erro = d_res_fd_jac.sum_abs_all();
    if (erro>1e-10) printf("error %g\n",erro);
#endif	

  } catch (GenericError e) {
    set_error(1);
    return;
  } 
  ADVABSO_CURRENT_ELEMENT = INT_MAX;
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  close();
} catch (GenericError e) {
  set_error(1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define ICONE(j,k) (icone[nel*(j)+(k)]) 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Load local node coordinates in local vector
void LagrangeMult::
get_xloc(FastMat2 &xloc,FastMat2 &Hloc) {
  const int 
    &nu = nodedata_m->nu,
    &ndim = nodedata_m->ndim;
  if (xloc.size()!=nel*ndim || Hloc.size()!=nel*nu) {
    xloc.resize(2,nel,ndim);
    Hloc.resize(2,nel,nu);
  }
  element.node_data(nodedata_m,xloc.storage_begin(),
		    Hloc.storage_begin());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Load local node coordinates in local vector
void LagrangeMult::
get_connect(vector<int> &node_list) {
  node_list.resize(nel);
  element_connect(element,&node_list[0]);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double LagrangeMult::
get_time() { return time_m; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LagrangeMult::
initialize() {
  lm_initialize();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LagrangeMult::
get_old_state(FastMat2 &Uold) {
  Uold.set(element.vector_values(*stateo));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LagrangeMult::
set_ldf(FastMat2 &ldf_user,
	vector<double> &ldf) { }
