//__INSERT_LICENSE__
//$Id: fem.cpp,v 1.15.10.1 2007/02/19 20:23:56 mstorti Exp $

#include <time.h>
#include <stdarg.h>

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/utils.h>
#include <src/getprop.h>
#include <src/pfmat.h>
#include <src/dvecpar.h>
#include <src/h5utils.h>
//#include <src/fastmat2.h>

MPI_Comm PETSCFEM_COMM_WORLD=0;
int MY_RANK=0, SIZE=1;

int PetscFemInitialize(int *argc,char ***args,
                       const char file[],const char help[]) {
  int PFUNUSED ierr;
  ierr = PetscInitialize(argc,args,file,help);
  PETSCFEM_COMM_WORLD = PETSC_COMM_WORLD;
  ierr = MPI_Comm_size(PETSCFEM_COMM_WORLD,&SIZE);
  ierr = MPI_Comm_rank(PETSCFEM_COMM_WORLD,&MY_RANK);

  // Just to resolve some linking problems.
  // The problem is that as the binaries are linked
  // statically, and the function `dvector_clone_parallel()'
  // is not called explicitly in PETSc-FEM, but it is
  // instantiated explicitly in `dvecpar.cpp'.
  // Then it seems that when the binary advdif_O.bin is linked
  // the instantiation are stripped off. Then when the function
  // is called from a dynamically linked file, it is not found
  // in the binary. This only happens with `advdif_O', not with
  // `ns', neither with the debugger version of `advdif'.
  // And only happens in some systems, for instance it happens in
  // aquiles (Fedora 9), but not in minerva (Fedora 7), neither
  // apparently on more recent versions of Fedora. 
  dvector<double> ppd;
  // It seems that now the compiler requires that the cloned array must not be empty
  ppd.a_resize(1,2);
  dvector_clone_parallel(ppd);
  dvector<int> ppi;
  ppi.a_resize(1,2);
  dvector_clone_parallel(ppi);
  dvector<float> ppf;
  ppf.a_resize(1,2);
  dvector_clone_parallel(ppf);
  dvector<char> ppc;
  ppc.a_resize(1,2);
  dvector_clone_parallel(ppc);
  h5_dvector_read(NULL,NULL,ppd);

  FastMat2::init();

  PETSCFEM_ASSERT0(ierr==0,"Initialization error");  
  return 0;
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "zeroe_mat"
int zeroe_mat(Mat A,int & ass_flag) {
  if (ass_flag==0) {
    ass_flag=1;
    return 0;
  } else {
    int ierr = MatZeroEntries(A);
    return ierr;
  }
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "compute_prof"
int compute_prof(Darray *da,Dofmap *dofmap,int myrank,
		 Mat *A, int debug_compute_prof=0) {
  int k,k1,k2,neqp,keq,leq,pos,sumd=0,sumdcorr=0,sumo=0,ierr;
  Node *nodep;
  k1=dofmap->startproc[myrank];
  neqp=dofmap->neqproc[myrank];
  k2=k1+neqp-1;

  // printf("on proc %d k1=%d k2=%d\n",myrank,k1,k2);
  int *d_nnz,*o_nnz,diag_ok;
  d_nnz = new int[neqp];
  o_nnz = new int[neqp];
  for (k=0;k<neqp;k++) {
    d_nnz[k]=0;
    o_nnz[k]=0;
    keq=k1+k;
    if (debug_compute_prof) printf("-------- keq = %d: ",keq);
    pos=keq;
    diag_ok=0;
    while (1) {
      nodep = (Node *)da_ref(da,pos);
      if (nodep->next==-1) break;
      leq = nodep->val;
      if (debug_compute_prof) printf("%d ",leq);
      if (leq==keq) diag_ok=1;
      if (k1<=leq && leq<=k2) {
	d_nnz[k]++;
      } else {
	o_nnz[k]++;
      }	
      pos = nodep->next;
    }
    if (debug_compute_prof) printf("     --    d_nnz %d   o_nnz %d\n",d_nnz[k],o_nnz[k]);
    //d_nnz[k] += 1;
    //o_nnz[k] += 1;
    // To add room for the diagonal entry, even if it doesn't exist
    sumd += d_nnz[k];
    if (!diag_ok) {
      //      printf("No diagonal element for dof %d in proc %d\n",keq,myrank);
      // fixme:= No hago correccion para elementos diagonales
      //d_nnz[k]+=1; 
    }
    sumdcorr += d_nnz[k];
    sumo += o_nnz[k];
  }

  // deallocate darray
  // wait_from_console("antes de destruir da ");
  da_destroy(da);
  // wait_from_console("despues de destruir da ");

  double avo,avd, PFUNUSED avdcorr;
  avo = double(sumo)/double(neqp);
  avd = double(sumd)/double(neqp);
  avdcorr = double(sumdcorr)/double(neqp);
  PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,
			  "On processor %d,\n"
			  "       diagonal block terms: %d, (%f av.)\n"
			  // Corrected does not make sense anymore
			  // since it seems that PETSc does not need
			  // the diagonal terms. 
			  // "                (corrected): %d, (%f av.)\n"
			  "   off diagonal block terms: %d, (%f av.)\n",
			  // myrank,sumd,avd,sumdcorr,avdcorr,sumo,avo);
			  myrank,sumd,avd,sumo,avo);
  PetscSynchronizedFlush(PETSCFEM_COMM_WORLD);
  
  // Create matrices
  int neq=dofmap->neq;
  ierr =  MatCreateMPIAIJ(PETSCFEM_COMM_WORLD,dofmap->neqproc[myrank],
			  dofmap->neqproc[myrank],neq,neq,
			  PETSC_NULL,d_nnz,PETSC_NULL,o_nnz,A); CHKERRA(ierr);
  delete[] d_nnz;
  delete[] o_nnz;
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void petscfem_printf(const char *,va_list )"
void petscfem_printf(const char *templ,va_list list) {
  int myrank;
  // va_list list;
  // va_start(list,templ);
  MPI_Comm_rank(PETSCFEM_COMM_WORLD,&myrank);
  if (myrank==0) vprintf(templ,list);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void petscfem_error(const char *,...)"
void petscfem_error(const char *templ,...) {
  va_list list;
  va_start(list,templ);
  petscfem_printf(templ,list);
  PetscFinalize();
  abort();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void petscfem_warn(const char *templ,...) {
  va_list list;
  va_start(list,templ);
  petscfem_printf(templ,list);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void petscfem_assert(int cond, const char *templ,...) {
  if (!cond) {
    va_list list;
    va_start(list,templ);
    petscfem_error(templ,list);
  }
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int opt_read_vector(Mesh *,Vec , Dofmap *,int )"
int opt_read_vector(Mesh *mesh,Vec x, Dofmap *dofmap,int myrank) {

  char *ini_name;
  int ierr;
  mesh->global_options->get_entry("initial_state",ini_name);
  if (ini_name!=NULL) {
    ierr = read_vector(ini_name,x,dofmap,myrank);
  } else {
    double scal = 0.;
    ierr = VecSet(x,scal); CHKERRA(ierr);
  }
  return ierr;
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int opt_read_vector(Mesh *,Vec , Dofmap *,int )"
int opt_read_vector(Mesh *mesh,Vec x, Dofmap *dofmap,int myrank) {
  string ini_name;
  int ierr;
  ini_name="";
  ierr = get_string(mesh->global_options,"initial_state",ini_name,1);
  if (ini_name!=string("")) {
    ierr = read_vector(ini_name.c_str(),x,dofmap,myrank);
  } else {
    double scal = 0.;
    ierr = VecSet(x,scal); CHKERRA(ierr);
  }
  return ierr;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void petscfem_check_par_err(int ierro,GenericError &ge,int myrank) {
  MPI_Bcast(&ierro,1,MPI_INT,0,PETSCFEM_COMM_WORLD);
  PETSCFEM_ASSERT(!ierro,"%s",ge.c_str());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void detj_error(double &detJaco,int elem) {
  printf("Jacobian of element %d is negative or null\n"
	 " Jacobian: %f\n",elem,detJaco);
  detJaco = -detJaco;
  if (detJaco==0.) detJaco = 1.0;
}

