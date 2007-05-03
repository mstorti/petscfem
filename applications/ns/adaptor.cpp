//__INSERT_LICENSE__
//$Id: adaptor.cpp,v 1.18.2.2 2007/03/20 17:41:23 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"

extern TextHashTable *GLOBAL_OPTIONS;
   
#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void adaptor::export_vals(adaptor::ArgHandle h,
                          double *vals,int s) {
  PETSCFEM_ASSERT0(elem>=0 && ielh>=0,
                   "Current element not set. "
                   "`export_vals()' was called probably from outside \n"
                   "adaptor' element loop\n");
  int index = h.index();
  arg_data &arg = (*arg_data_vp)[index];
  PETSCFEM_ASSERT0(!(h==NullArgHandle),
                   "Invalid handle");
  PETSCFEM_ASSERT(index>=0 && index<int(nargs()),
                  "Invalid handle index, index %d, nargs\n",
                  index,nargs());  
  double *retval = arg.retval;
  int rowsz=-1;
  if (arg.options & OUT_VECTOR) rowsz=nel*ndof;
  else if (arg.options &OUT_MATRIX) 
    rowsz=nel*ndof*nel*ndof;
  else {
    PETSCFEM_ERROR("Invalid argument for output. \n"
                   "arginfo \"%s\", index %d\n",
                   arg.arginfo.c_str(),index);  
  }
  if (s>=0) {
    PETSCFEM_ASSERT(s==rowsz,
                    "Not exporting correct size. \n"
                    "exported %d, required %d\n",s,rowsz);  
  }
  memcpy(retval+ielh*rowsz,vals,rowsz*sizeof(double));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void adaptor::export_vals(ArgHandle h,FastMat2 &a) {
  export_vals(h,a.storage_begin(),a.size());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
adaptor::ArgHandle adaptor::NullArgHandle;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
size_t adaptor::nargs() const {
  return arg_data_vp->size();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
bool adaptor::ArgHandle::
operator==(ArgHandle b) const {
  return index_m==b.index_m;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
adaptor::ArgHandle 
adaptor::get_arg_handle(const string &key,
                        const char* errmess) const {
  ArgHandle handle = NullArgHandle;
  if (arg_data_vp) {
    for (unsigned int j=0; j<arg_data_vp->size(); j++) {
      if ((*arg_data_vp)[j].arginfo == key) {
        handle = ArgHandle(j);
        break;
      }
    }
  }
  if (errmess && handle==NullArgHandle)
    petscfem_error(errmess);
  return handle;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
adaptor::adaptor() : elem_init_flag(0), 
                     use_fastmat2_cache(1),
                     arg_data_vp(NULL),
                     elem(-1)
                     { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void adaptor::after_assemble(const char *jobinfo) {
  GET_JOBINFO_FLAG(comp_mat_res);
  if (comp_mat_res && !elem_init_flag) elem_init_flag=1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// modif nsi_tet
#undef __FUNC__
#define __FUNC__ "nsi_tet_les_fm2::assemble"
int adaptor::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
		      Dofmap *dofmap,const char *jobinfo,int myrank,
		      int el_start,int el_last,int iter_mode,
		      const TimeData *time_) {

  int kloc,node;
  arg_data_vp = &arg_data_v;

  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_mat_res);
  GET_JOBINFO_FLAG(comp_res);

  assert(!comp_res);

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retvalmat,iele,j,nel,k,ndof,p,nel,q,ndof)

  int ierr=0;
  // PetscPrintf(PETSC_COMM_WORLD,"entrando a nsilesther\n");

#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define ELEMIPROPS_ADD(j,k) VEC2(elemiprops_add,j,k,neliprops_add)

  //o Number of Gauss points.
  TGETOPTDEF_ND(thash,int,npg,0);
  // ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  TGETOPTDEF_ND(thash,int,ndim,0); //nd
  //o Use #arg-handles# for manipulation of argumentes
  //  from/to `adaptor'. 
  TGETOPTDEF(thash,int,use_arg_handles,0);
  //o Use caches for FastMat2 matrices
  TGETOPTDEF(thash,int,use_fastmat2_cache,1);

  PETSCFEM_ASSERT(npg>=0,"npg should be non-negative, npg %d\n",npg);  
  PETSCFEM_ASSERT(ndim>=0,"ndim should be non-negative, ndim %d\n",ndim);  

  TGETOPTDEF_ND(thash,int,ndimel,ndim);
  PETSCFEM_ASSERT(ndimel>=0 && ndimel<=ndim,
                  "Incorrect value for `ndimel': \n",ndimel);
  int nen = nel*ndof;

  // Unpack nodedata
  int nu=nodedata->nu;

  // Get arguments from arg_list
  double *locst=NULL,*locst2=NULL,*retval=NULL,*retvalmat=NULL;

  if (!use_arg_handles) {
    if (comp_mat) 
      retvalmat = arg_data_v[0].retval;
    if (comp_mat_res) {
      int ja=0;
      locst = arg_data_v[ja++].locst;
      locst2 = arg_data_v[ja++].locst;
      retval = arg_data_v[ja++].retval;
      retvalmat = arg_data_v[ja++].retval;
      ja++;
      glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
      alpha = glob_param->alpha;
    } 
  } else {

#define GET_INDEX(key)                                  \
    handle = get_arg_handle(key,                        \
          "can't get arg handle for key " key "\n");    \
    index = handle.index();

    ArgHandle handle;
    int index;
    if (comp_mat) {
      GET_INDEX("A"); retvalmat = arg_data_v[index].retval;
    } else if (comp_mat_res) {
      GET_INDEX("state"); locst = arg_data_v[index].locst;
      GET_INDEX("state_old"); locst2 = arg_data_v[index].locst;
      GET_INDEX("res"); retval = arg_data_v[index].retval;
      GET_INDEX("A"); retvalmat = arg_data_v[index].retval;
      GET_INDEX("glob_param"); 
      glob_param = (GlobParam *)(arg_data_v[index].user_data);
      alpha = glob_param->alpha;
    } 
  }

  if (comp_mat_res) {
    rec_Dt = 1./glob_param->Dt;
    if (glob_param->steady) rec_Dt=0.;
  }

  // allocate local vecs
  nen = nel*ndof;
  FastMat2 veccontr(2,nel,ndof),veccontrp(2,nel,ndof),
    xloc(2,nel,ndim),
    locstate(2,nel,ndof), locstatep(2,nel,ndof), 
    locstate2(2,nel,ndof), 
    matlocf(4,nel,ndof,nel,ndof),
    matlocf_fdj(4,nel,ndof,nel,ndof),
    matloc_prof(4,nel,ndof,nel,ndof), tmp;
    
  nH = nu-ndim;
  Hloc.resize(2,nel,nH);

  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  //GPdata gp_data(geom,ndim,nel,npg);
  GPdata gp_data(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);
  //o Compute a Finite Difference Jacobian (FDJ) for each element. 
  TGETOPTDEF(thash,int,jacobian_fdj_compute,0);
  //o Scale of perturbation for computation of FDJ's. 
  TGETOPTDEF(thash,double,jacobian_fdj_epsilon,1e-3);
  //o Print computed FDJ for comparison with analytic Jacobian. 
  TGETOPTDEF(thash,int,jacobian_fdj_print,1);
  //o Use FDJ for computations. 
  TGETOPTDEF(thash,int,use_jacobian_fdj,0);
  //o Compute variables #H#
  TGETOPTDEF(thash,int,compute_H_fields,0);

  if (use_jacobian_fdj) 
    jacobian_fdj_compute = 1;
  
#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

  // Memory allocation and initialization 
  shape.resize(2,nel,npg);
  dshapexi.resize(3,ndimel,nel,npg);
  wpg.resize(1,npg);

  for (int ipg=0; ipg<npg; ipg++) {
    wpg.setel(WPG,ipg+1);
    shape.ir(2,ipg+1).set(SHAPE);
    dshapexi.ir(3,ipg+1).set(DSHAPEXI);
  }
  shape.rs();
  dshapexi.rs();

  matloc_prof.set(1.);

  if (comp_mat_res && !elem_init_flag) {
    for (elem=el_start; elem<=el_last; elem++) {
      if (!compute_this_elem(elem,this,myrank,iter_mode)) continue;
      element_init();
    }
  }
  if (error_code()) return error_code();

  // Users may use `init()' in order to perform calculations
  // *outside* the element loop
  if (comp_mat_res) init();

  before_chunk(jobinfo);
  FastMatCacheList cache_list;
  if (use_fastmat2_cache) 
    FastMat2::activate_cache(&cache_list);

  ielh=-1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();
    ielh++; elem=k;
    // load_props(propel,elprpsindx,nprops,&(ELEMPROPS(k,0)));

    // Load local node coordinates in local vector
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(k,kloc);
      xloc.ir(1,kloc+1).set(&NODEDATA(node-1,0));
      if (nH>0 && compute_H_fields)
	Hloc.ir(1,kloc+1).set(&NODEDATA(node-1,0)+ndim);
    }
    Hloc.rs();
    xloc.rs();

    if (comp_mat_res) {
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));
    }

    matlocf.set(0.);
    veccontr.set(0.);

    if(comp_mat && !use_arg_handles)
      matloc_prof.export_vals(&(RETVALMAT(ielh,0,0,0,0)));

    if (comp_mat_res) {
      // Users have to implement this function with the physics of
      // the problem.
      element_connector(xloc,locstate2,locstate,veccontr,matlocf);
      if (!use_arg_handles)
        veccontr.export_vals(&(RETVAL(ielh,0,0)));

      if (jacobian_fdj_compute) {
	double epsil = jacobian_fdj_epsilon;
	matlocf_fdj.set(0.).reshape(2,nen,nen);
	veccontr.reshape(1,nen);
	locstate.reshape(1,nen);
	double afact = -1./(alpha*epsil);
	for (int j=1; j<=nen; j++) {
	  locstatep.reshape(1,nen).set(locstate)
	    .addel(epsil,j).reshape(2,nel,ndof);
	  veccontrp.reshape(2,nel,ndof);
	  element_connector(xloc,locstate2,locstatep,veccontrp,matlocf);
	  veccontrp.reshape(1,nen);
	  matlocf_fdj.ir(2,j).set(veccontrp)
	    .rest(veccontr).scale(afact).rs();
	}
	veccontr.reshape(2,nel,ndof);
	locstate.reshape(2,nel,ndof);
	if (jacobian_fdj_print) {
	  matlocf.reshape(4,nel,ndof,nel,ndof);
	  tmp.ctr(matlocf,2,1,4,3);
	  tmp.print(nen,"analytic jacobian: ");
	  matlocf_fdj.reshape(4,nel,ndof,nel,ndof);
	  tmp.ctr(matlocf_fdj,2,1,4,3);
	  tmp.print(nen,"FD jacobian: ");
	  matlocf.rs();
	  matlocf_fdj.rs();
	}
	if (use_jacobian_fdj) 
	  matlocf_fdj.reshape(4,nel,ndof,nel,ndof);
	  matlocf.set(matlocf_fdj);
	  matlocf_fdj.rs();
      }
      if (!use_arg_handles)
        matlocf.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }
  }
  ielh = elem = -1;
  after_chunk(jobinfo);
  FastMat2::void_cache();
  FastMat2::deactivate_cache();

  if (comp_mat_res) clean();
  arg_data_vp = NULL;
  return 0;
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
