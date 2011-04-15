//__INSERT_LICENSE__
//$Id merge-with-petsc-233-55-g52bd457 Fri Oct 26 13:57:07 2007 -0300$

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
void adaptor::get_vals(adaptor::ArgHandle h,
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
  double *locst = arg.locst;
  int rowsz=-1;
  PETSCFEM_ASSERT(arg.options & IN_VECTOR,
                  "Arg handle (index=%d,key=\"%s\") should be an IN_VECTOR",
                  index,arg.arginfo.c_str());
  rowsz=nel*ndof;
  if (s>=0) {
    PETSCFEM_ASSERT(s==rowsz,
                    "Not given correct size. \n"
                    "got %d, required %d\n",s,rowsz);  
  }
  memcpy(vals,locst+ielh*rowsz,rowsz*sizeof(double));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void adaptor::get_vals(ArgHandle h,FastMat2 &a) {
  get_vals(h,a.storage_begin(),a.size());
}

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
  if (arg.options & OUT_VECTOR) 
    rowsz=nel*ndof;
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
                     jpert(-1), 
                     elem(-1)
                     { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void adaptor::after_assemble(const char *jobinfo) {
  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_mat_res);
  if (!elem_init_flag &&
      (comp_res || comp_mat || comp_mat_res)
      ) elem_init_flag=1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int adaptor::prtb_index() { return jpert; }

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

  GET_JOBINFO_FLAG(comp_prof);
  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_mat_res);
  if (comp_mat_res) comp_res = comp_mat = 1;
  
#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retvalmat,iele,j,nel,k,ndof,p,nel,q,ndof)

  int ierr=0;
  // PetscPrintf(PETSCFEM_COMM_WORLD,"entrando a nsilesther\n");

#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define ELEMIPROPS_ADD(j,k) VEC2(elemiprops_add,j,k,neliprops_add)

  //o Number of Dimensions.
  TGETOPTDEF_ND(thash,int,ndim,nodedata->ndim);
  PETSCFEM_ASSERT(ndim>=0,
		  "ndim should be non-negative, ndim %d\n",ndim);  
  TGETOPTDEF_ND(thash,int,ndimel,ndim);
  PETSCFEM_ASSERT(ndimel>=0 && ndimel<=ndim,
                  "Incorrect value for `ndimel': \n",ndimel);

  //o Use #arg-handles# for manipulation of argumentes
  //  from/to `adaptor'. 
  TGETOPTDEF(thash,int,use_arg_handles,0);
  //o Use caches for FastMat2 matrices
  TGETOPTDEF(thash,int,use_fastmat2_cache,1);

  int nen = nel*ndof;

  // Unpack nodedata
  int nu=nodedata->nu;

  // Get arguments from arg_list
  double *locst=NULL,*locst2=NULL,*retval=NULL,*retvalmat=NULL;

  if (!use_arg_handles) {
    if (comp_prof) 
      retvalmat = arg_data_v[0].retval;
    else if (comp_res || comp_mat) {
      int ja=0;
      locst = arg_data_v[ja++].locst;
      locst2 = arg_data_v[ja++].locst;
      retval = arg_data_v[ja++].retval;
      if (comp_mat) {
	retvalmat = arg_data_v[ja++].retval;
      }
      ja++;
      glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
      alpha = glob_param->alpha;
      if (glob_param->steady) rec_Dt=0.;
      else rec_Dt = 1./glob_param->Dt;
    } 
  } else {

#define GET_INDEX(key)                                  \
    handle = get_arg_handle(key,                        \
          "can't get arg handle for key " key "\n");    \
    index = handle.index();

    ArgHandle handle;
    int index;
    if (comp_prof) {
      GET_INDEX("A"); retvalmat = arg_data_v[index].retval;
    } else if (comp_res || comp_mat) {
      GET_INDEX("state"); locst = arg_data_v[index].locst;
      GET_INDEX("state_old"); locst2 = arg_data_v[index].locst;
      GET_INDEX("res"); retval = arg_data_v[index].retval;
      if (comp_mat) {
	GET_INDEX("A"); retvalmat = arg_data_v[index].retval;
      }
      GET_INDEX("glob_param"); 
      glob_param = (GlobParam *)(arg_data_v[index].user_data);
      alpha = glob_param->alpha;
      if (glob_param->steady) rec_Dt=0.;
      else rec_Dt = 1./glob_param->Dt;
    } 
  }
  // allocate local vecs
  nen = nel*ndof;
  FastMat2 
    veccontr(2,nel,ndof),
    veccontra(2,nel,ndof),
    veccontrp(2,nel,ndof),
    veccontrm(2,nel,ndof),
    xloc(2,nel,ndim),
    locstate(2,nel,ndof), locstatep(2,nel,ndof), 
    locstate2(2,nel,ndof), 
    matlocf(4,nel,ndof,nel,ndof),
    matlocfa(4,nel,ndof,nel,ndof),
    matlocfp(4,nel,ndof,nel,ndof),
    matlocf_fdj(4,nel,ndof,nel,ndof),
    matloc_prof(4,nel,ndof,nel,ndof), tmp;
    
  nH = nu-ndim;
  Hloc.resize(2,nel,nH);

  //o Compute a Finite Difference Jacobian (FDJ) for each element. 
  TGETOPTDEF(thash,int,jacobian_fdj_compute,0);
  //o Scale of perturbation for computation of FDJ's. 
  TGETOPTDEF(thash,double,jacobian_fdj_epsilon,1e-7);
  //o Print computed FDJ for comparison with analytic Jacobian. 
  TGETOPTDEF(thash,int,jacobian_fdj_print,0);
  //o Use FDJ for computations. 
  TGETOPTDEF(thash,int,use_jacobian_fdj,0);
  //o Uses second order formular for the evaluation
  // of numerical Jacobians (bur requires twice the number of
  // evaluations of residuals). In the low order version the
  // Jacobian is approximated by #J_k = (R(x+eps*e_k)-R(x))/eps#,
  // whereas the higher order version is
  // #J_k = (R(x+eps*e_k)-R(x-eps*e_k))/(2*eps)#. #ek# is
  // the unit vector along the #k# axis. 
  TGETOPTDEF(thash,int,use_high_prec,0);
  //o Compute variables #H#
  TGETOPTDEF(thash,int,compute_H_fields,0);

  if (use_jacobian_fdj || jacobian_fdj_print) 
    jacobian_fdj_compute = 1;
  
#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

  GPdata gp_data;
  //o Use #shape()#, #shapexi()# (FEM interpolation)
  TGETOPTDEF(thash,int,fem_interpolation,1);
  if (fem_interpolation) {
    //o Type of element geometry to define Gauss Point data
    TGETOPTDEF_S(thash,string,geometry,default);
    //o Number of Gauss points.
    TGETOPTDEF_ND(thash,int,npg,-1);
    if (geometry == "default") GPdata::get_default_geom(ndimel, nel, geometry);
    if (npg      == -1       ) GPdata::get_default_npg(geometry, npg);
    PETSCFEM_ASSERT0(geometry!="default","could not determine geometry");
    PETSCFEM_ASSERT(npg>=0,"npg should be non-negative, npg %d\n",npg);

    gp_data.init(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);
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
  }

  if ((comp_res || comp_mat) && !elem_init_flag) {
    for (elem=el_start; elem<=el_last; elem++) {
      if (!compute_this_elem(elem,this,myrank,iter_mode)) continue;
      element_init();
    }
  }
  if (error_code()) return error_code();

  // Users may use `init()' in order to perform calculations
  // *outside* the element loop
  if (comp_res || comp_mat) init();

  before_chunk(jobinfo);
  FastMatCacheList cache_list;
  if (use_fastmat2_cache) 
    FastMat2::activate_cache(&cache_list);

  EVAL_RES = EVAL_MAT = false;
  if (comp_res) EVAL_RES = true;
  if (comp_mat) EVAL_MAT = true;
  
  if (comp_mat && jacobian_fdj_compute) EVAL_RES = true;
  if (comp_mat && jacobian_fdj_print)   EVAL_MAT = true;


  ielh=-1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();
    ielh++; elem=k;
    // load_props(propel,elprpsindx,nprops,&(ELEMPROPS(k,0)));

    if(comp_prof && !use_arg_handles) {
      matloc_prof.set(1.);
      matloc_prof.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }
    
    // Load local node coordinates in local vector
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(k,kloc);
      xloc.ir(1,kloc+1).set(&NODEDATA(node-1,0));
      if (nH>0 && compute_H_fields)
	Hloc.ir(1,kloc+1).set(&NODEDATA(node-1,0)+ndim);
    }
    Hloc.rs();
    xloc.rs();

    if (comp_res || comp_mat) {
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));
    }

    if (comp_res || comp_mat) {

      if (EVAL_RES) veccontr.set(0.);
      if (EVAL_MAT) matlocf.set(0.);

      // Users have to implement this function 
      // with the physics of the problem.
      jpert=0;
      element_connector(xloc,locstate2,locstate,veccontr,matlocf);
      veccontra.set(0.0);
      matlocfa.set(0.0);
      jpert=-1;
      element_connector_analytic(xloc,locstate2,locstate,veccontra,matlocfa);

      if (comp_res) {
        veccontra.add(veccontr);
        veccontra.export_vals(&(RETVAL(ielh,0,0)));
      }

      if (comp_mat && jacobian_fdj_compute) {
	double epsil = jacobian_fdj_epsilon;
	double afact = -1./(alpha*epsil);
	matlocf_fdj.set(0.).reshape(2,nen,nen);
	veccontr.reshape(1,nen);
	locstate.reshape(1,nen);
	for (jpert=1; jpert<=nen; jpert++) {
	  locstatep.reshape(1,nen).set(locstate)
	    .addel(epsil,jpert).reshape(2,nel,ndof);
	  veccontrp.reshape(2,nel,ndof).set(0.0);
	  element_connector(xloc,locstate2,locstatep,
                            veccontrp,matlocfp);
	  veccontrp.reshape(1,nen);
          if (!use_high_prec) {
            matlocf_fdj.ir(2,jpert).set(veccontrp)
              .minus(veccontr).scale(afact).rs();
          } else {
            locstatep.reshape(1,nen).set(locstate)
              .addel(-epsil,jpert).reshape(2,nel,ndof);
            veccontrm.reshape(2,nel,ndof).set(0.0);
            element_connector(xloc,locstate2,locstatep,
                              veccontrm,matlocfp);
            veccontrm.reshape(1,nen);
            double af = -1.0/(2.0*alpha*epsil);
            matlocf_fdj.ir(2,jpert).set(veccontrp)
              .minus(veccontrm).scale(af).rs();
          }
	}
        jpert=-1;
	matlocf_fdj.reshape(4,nel,ndof,nel,ndof);
	veccontr.reshape(2,nel,ndof);
	locstate.reshape(2,nel,ndof);
	if (jacobian_fdj_print) {
	  tmp.ctr(matlocf,2,1,4,3);
	  tmp.print(nen,"analytic jacobian: ");
	  tmp.ctr(matlocf_fdj,2,1,4,3);
	  tmp.print(nen,"FD jacobian: ");
	}
	if (use_jacobian_fdj)
	  matlocf.set(matlocf_fdj);
      }
      
      // if (comp_mat && !use_arg_handles) {
      if (comp_mat) {
        matlocfa.add(matlocf);
        matlocfa.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
      }
    }
  }

  ielh = elem = -1;
  after_chunk(jobinfo);
  FastMat2::void_cache();
  FastMat2::deactivate_cache();

  if (comp_res || comp_mat) clean();
  arg_data_vp = NULL;
  return 0;
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
