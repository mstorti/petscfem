//__INSERT_LICENSE__
/* $Id: nssupr.cpp,v 1.12 2003/07/06 15:10:18 mstorti Exp $ */

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <applications/ns/nsi_tet.h>
#include <applications/ns/nssup.h>

extern TextHashTable *GLOBAL_OPTIONS;

#ifdef ROSI_COUPLING_MODULE
#warning Compiling with ROSI_COUPLING_MODULE enabled
double AXIAL_ACCELERATION=DBL_MAX, GLOBAL_TIME;
extern int MY_RANK,SIZE;
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void ns_sup_res::lm_initialize() {
#ifdef ROSI_COUPLING_MODULE
  int ierr;
  //o Is being PETSc-FEM called from ROSI?
  TGETOPTDEF_ND(thash,int,called_from_rosi,0);
  // This is tricky :-)
  // `AXIAL_ACCELERATION' may have been set by an elemset like
  // `nsi_tet_keps_rot' in some of the processors, but may be
  // not in ALL processors!!
  // The rank where the `AXIAL_ACCELERATION' has been set. If it was not
  // set, we set it to -1. 
  int axial_accel_rank_a = (AXIAL_ACCELERATION!=DBL_MAX ? MY_RANK : -1);
#if 0
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] axial_accel_rank_a: %d\n",
			  MY_RANK,axial_accel_rank_a);
  PetscSynchronizedFlush(PETSC_COMM_WORLD); 
#endif
  int axial_accel_rank;
  // The we make an ALlreduce and get one of the processors where
  // the acceleration has been set.
  MPI_Allreduce(&axial_accel_rank_a,&axial_accel_rank,1,MPI_INT,
		MPI_MAX,PETSC_COMM_WORLD);
  if (axial_accel_rank>=0)
    // And we make a broadcast from that processor
    MPI_Bcast(&AXIAL_ACCELERATION,1,MPI_DOUBLE,axial_accel_rank,PETSC_COMM_WORLD);
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void ns_sup_res::init() {
  int ierr=0;
  //o Gravity acceleration. 
  TGETOPTDEF_ND(thash,double,gravity,1.);
  //o Gravity acceleration. 
  TGETOPTDEF_ND(thash,double,rho,1.);
  //o Dimension of problem
  TGETOPTDEF_ND(thash,int,ndim,0);
  assert(ndim>0);
  p_indx = ndim+1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
void ns_sup_res::res(int k,FastMat2 & U,FastMat2 & r,
		       FastMat2 & w,FastMat2 & jac) {
  double p, eta;

  // We have (for 3D) first node: {u,v,w,p}
  //                  second node: {eta,lambda,*,*}
  assert(nel==2);
  // Discard k and eps eqs.
  w.set(0.).setel(1.,1,ndim+1,1);
  
#ifdef ROSI_COUPLING_MODULE
  double total_axial_acc = gravity;
  // Verify that AXIAL_ACCELERATION has been set (probably in `rosi_hook')
  if (called_from_rosi) {
    assert(AXIAL_ACCELERATION!=DBL_MAX);
    total_axial_acc += AXIAL_ACCELERATION;
    // We could have problems with the free surface if the axal acceleration is too low
    if (total_axial_acc < 0.3 * gravity) total_axial_acc = 0.3 * gravity;
  }
#define gravity total_axial_acc
#endif

  p = U.get(1,p_indx);
  eta = U.get(2,1);
  //#define DEBUG_AXIAL_ACC
#if DEBUG_AXIAL_ACC
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"k %d, gravity %f, rho %f, p %f, eta %f\n",
			  k, gravity, rho, p, eta);
#endif
  r.setel(p-rho*gravity*eta,1);
  jac.setel(1.,1,1,p_indx).setel(-rho*gravity,1,2,1);
}

void ns_sup_res::close() {
#if DEBUG_AXIAL_ACC
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
#endif
}
