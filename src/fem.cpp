//__INSERT_LICENSE__
//$Id: fem.cpp,v 1.11 2003/07/03 04:32:11 mstorti Exp $

#include <time.h>
#include <stdarg.h>

#include "fem.h"
#include "readmesh.h"
#include "utils.h"
#include "getprop.h"
#include "pfmat.h"

int MY_RANK, SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "zeroe_mat"
int zeroe_mat(Mat A,int & ass_flag) {
  if (ass_flag==0) {
    ass_flag=1;
    return 0;
  } else {
    double scal = 0;
    int ierr = MatScale(&scal,A);
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

  double avo,avd,avdcorr;
  avo = double(sumo)/double(neqp);
  avd = double(sumd)/double(neqp);
  avdcorr = double(sumdcorr)/double(neqp);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			  "On processor %d,\n"
			  "       diagonal block terms: %d, (%f av.)\n"
			  // Corrected does not make sense anymore
			  // since it seems that PETSc does not need
			  // the diagonal terms. 
			  // "                (corrected): %d, (%f av.)\n"
			  "   off diagonal block terms: %d, (%f av.)\n",
			  // myrank,sumd,avd,sumdcorr,avdcorr,sumo,avo);
			  myrank,sumd,avd,sumo,avo);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  
  // Create matrices
  int neq=dofmap->neq;
  ierr =  MatCreateMPIAIJ(PETSC_COMM_WORLD,dofmap->neqproc[myrank],
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
//    va_list list;
//    va_start(list,templ);
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
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
  exit(0);
}

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
    ierr = VecSet(&scal,x); CHKERRA(ierr);
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
    ierr = VecSet(&scal,x); CHKERRA(ierr);
  }
  return ierr;
}
