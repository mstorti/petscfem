//__INSERT_LICENSE__
//$Id: ns-fracstep.cpp,v 1.9 2002/09/05 20:10:17 mstorti Exp $
 
#include <time.h>
#include <malloc.h>

#include "../../src/fem.h"
#include "../../src/readmesh.h"
#include "../../src/getprop.h"
#include "fracstep.h"
#include "nsi_tet.h"

static char help[] = "Basic finite element program.\n\n";

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

int MyKSPMonitor(KSP ksp,int n,double rnorm,void *dummy)
{
  Vec      x;
  int      ierr;

  /* 
     Build the solution vector
  */
  ierr = KSPBuildSolution(ksp,PETSC_NULL,&x); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,
	      "iteration %d KSP Residual norm %14.12e \n",n,rnorm);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "bless_elemset"
void bless_elemset(char *type,Elemset *& elemset) {
  //  SET_ELEMSET_TYPE(internal)
    SET_ELEMSET_TYPE(fracstep)
    SET_ELEMSET_TYPE(nsi_tet)
	{
	printf("not known elemset \"type\": %s\n",type);
	exit(1);
	}
}

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
#undef __FUNC__
#define __FUNC__ "main"
int main(int argc,char **args) {

// #define NS   // Navier Stokes eqs. by fractional step
#define NSI_TET   // Navier Stokes eqs. by SUPG_PSPG (Tezduyar et.al)

#ifdef NS
  Vec     x, dx, xold, dx_step, res;      /* approx solution, RHS, residual*/
  Mat     A_mom, A_poi, A_prj;            /* linear system matrix */
  SLES    sles_mom,sles_poi,sles_prj;     /* linear solver context */
  PC      pc_mom,pc_poi,pc_prj;           /* preconditioner context */
  KSP     ksp_mom,ksp_poi,ksp_prj;        /* Krylov subspace method context */
  double  norm, *sol, scal; /* norm of solution error */
  int     ierr, i, n = 10, col[3], its, flg, size, node,
    jdof, k, kk, nfixa,
    kdof, ldof, lloc, ndim, nel, nen, neq, nu,
    myrank;
  // nu:= dimension of the state vector per node
  PetscScalar  neg_one = -1.0, one = 1.0, value[3];
  PetscScalar *px;
  char fcase[FLEN+1];
  Darray *da;
  //Elemset *elemset;
  Dofmap *dofmap;
  dofmap = new Dofmap;
  Mesh *mesh;

  // elemsetlist =  da_create(sizeof(Elemset *));
  PetscInitialize(&argc,&args,(char *)0,help);

  // Get MPI info
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);

  //  if (size != 1) SETERRA(1,0,"This is a uniprocessor example only!");
  ierr = PetscOptionsGetString(PETSC_NULL,"-case",fcase,FLEN,&flg); CHKERRA(ierr);

  // Read data

//    FILE *fid;
//    char line[LINESIZE];
//    fid=fopen(fcase,"r");
//    if (! fid) {
//      PetscPrintf(PETSC_COMM_WORLD,"Couldn't open file %s!!\n",fcase);
//      exit(1);
//    }

  // Read the mesh
  //read_mesh(fid,&nodedata,dofmap,elemsetlist,neq,size,myrank,x,xseq);
  read_mesh(mesh,fcase,dofmap,neq,size,myrank,x,xseq);

  GETOPTDEF(double,tol_mom,5e-4);
  GETOPTDEF(double,tol_poi,5e-4);
  GETOPTDEF(double,tol_prj,5e-4);

  GETOPTDEF(int,nsave,10);
  GETOPTDEF(int,nstep,10000);
  GETOPTDEF(int,nnwt,10);

  // initialize vectors
  ierr = VecDuplicate(x,&xp); CHKERRA(ierr);
  ierr = VecDuplicate(x,&xold); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx_step); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx); CHKERRA(ierr);
  ierr = VecDuplicate(x,&res); CHKERRA(ierr);
  ierr = VecDuplicate(xseq,&xpseq); CHKERRA(ierr);
  ierr = VecDuplicate(xseq,&xoldseq); CHKERRA(ierr);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // COMPUTE ACTVE PROFILE ON A LIBRETTO DARRAY'S
  // initialize state vectors
  scal=0;
  ierr = VecSet(&scal,x); CHKERRA(ierr);
  ierr = VecSet(&scal,xp); CHKERRA(ierr);

  ierr = assemble(mesh,res,A_mom,da,xp,xpseq,x,xseq,
  		  dofmap,COMP_MAT_PROF,"comp_mat_mom_prof"); CHKERRA(ierr);
  // ierr = assemble(mesh,res,A_mom,da,xp,xpseq,x,xseq,
  // dofmap,COMP_MAT_PROF,"comp_mat_mom"); CHKERRA(ierr);

  // compute d_nnz and allocate PETSc matrix
  //debug_compute_prof=1;
  PetscPrintf(PETSC_COMM_WORLD,"Profile for momentum substep ---------\n");
  ierr = compute_prof(da,dofmap,myrank,&A_mom);
 
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // ierr = assemble(mesh,res,A_poi,da,xp,xpseq,x,xseq,
  // dofmap,COMP_FDJ_PROF,"comp_res_poi"); CHKERRA(ierr);
  ierr = assemble(mesh,res,A_poi,da,xp,xpseq,x,xseq,
		  dofmap,COMP_MAT_PROF,"comp_mat_poi_prof"); CHKERRA(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Profile for poisson substep ---------\n");
  ierr = compute_prof(da,dofmap,myrank,&A_poi);
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  ierr = assemble(mesh,res,A_prj,da,xp,xpseq,x,xseq,
		  dofmap,COMP_MAT_PROF,"comp_mat_prj"); CHKERRA(ierr);
  
  PetscPrintf(PETSC_COMM_WORLD,"Profile for projection substep ---------\n");
  ierr = compute_prof(da,dofmap,myrank,&A_prj);

  // Flags whether the matrices have been assembled or not
  int a_mom_assem=0,a_poi_assem=0,a_prj_assem=0;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // SLES para el paso predictor
  ierr = SLESCreate(PETSC_COMM_WORLD,&sles_mom); CHKERRA(ierr);
  ierr = SLESSetOperators(sles_mom,A_mom,
			  A_mom,DIFFERENT_NONZERO_PATTERN); CHKERRA(ierr);
  ierr = SLESGetKSP(sles_mom,&ksp_mom); CHKERRA(ierr);
  ierr = SLESGetPC(sles_mom,&pc_mom); CHKERRA(ierr);

  ierr = KSPSetType(ksp_mom,KSPGMRES); CHKERRA(ierr);
  ierr = KSPGMRESSetRestart(ksp_mom,50); CHKERRA(ierr);
  ierr = PCSetType(pc_mom,PCJACOBI); CHKERRA(ierr);
  ierr = KSPSetTolerances(ksp_mom,tol_mom,PETSC_DEFAULT,PETSC_DEFAULT,
         PETSC_DEFAULT); CHKERRA(ierr);
  ierr = KSPSetMonitor(ksp_mom,MyKSPMonitor,PETSC_NULL);

  // SLES para el paso de poisson
  //ierr = MatDuplicate(A_mom,MAT_COPY_VALUES,&A_poi); CHKERRA(ierr);
  ierr = SLESCreate(PETSC_COMM_WORLD,&sles_poi); CHKERRA(ierr);
  ierr = SLESSetOperators(sles_poi,A_poi,A_poi,
			  DIFFERENT_NONZERO_PATTERN); CHKERRA(ierr);
  ierr = SLESGetKSP(sles_poi,&ksp_poi); CHKERRA(ierr);
  ierr = SLESGetPC(sles_poi,&pc_poi); CHKERRA(ierr);

  ierr = KSPSetType(ksp_poi,KSPCG); CHKERRA(ierr);
  //ierr = KSPSetType(ksp_poi,KSPBICG); CHKERRA(ierr);
  //ierr = PCSetType(pc_poi,PCNONE); CHKERRA(ierr);
  ierr = PCSetType(pc_poi,PCJACOBI); CHKERRA(ierr);
  ierr = KSPSetTolerances(ksp_poi,tol_poi,PETSC_DEFAULT,PETSC_DEFAULT,
         PETSC_DEFAULT); CHKERRA(ierr);
  ierr = KSPSetMonitor(ksp_poi,MyKSPMonitor,PETSC_NULL);

  // SLES para el paso de projection
  //ierr = MatDuplicate(A_mom,MAT_COPY_VALUES,&A_prj); CHKERRA(ierr);
  ierr = SLESCreate(PETSC_COMM_WORLD,&sles_prj); CHKERRA(ierr);
  ierr = SLESSetOperators(sles_prj,A_prj,A_prj,
			  DIFFERENT_NONZERO_PATTERN); CHKERRA(ierr);
  ierr = SLESGetKSP(sles_prj,&ksp_prj); CHKERRA(ierr);
  ierr = SLESGetPC(sles_prj,&pc_prj); CHKERRA(ierr);

  ierr = KSPSetType(ksp_prj,KSPCG); CHKERRA(ierr);
  ierr = PCSetType(pc_prj,PCJACOBI); CHKERRA(ierr);
  ierr = KSPSetTolerances(ksp_prj,tol_prj,PETSC_DEFAULT,PETSC_DEFAULT,
         PETSC_DEFAULT); CHKERRA(ierr);
  ierr = KSPSetMonitor(ksp_prj,MyKSPMonitor,PETSC_NULL);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 1
      PetscViewer matlab;
      ierr = ViewerASCIIOpen(PETSC_COMM_WORLD,
			     "matns.m",&matlab); CHKERRA(ierr);
      ierr = ViewerSetFormat(matlab,
      			     PETSC_VIEWER_ASCII_MATLAB,"res");
      CHKERRA(ierr);
#endif 
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  char *ini_name;
  mesh->global_options->get_entry("initial_state",ini_name);
  if (ini_name!=NULL) {
    read_vector(ini_name,x,xseq,dofmap,myrank);
  } else {
    scal = 0.;
    ierr = VecSet(&scal,x); CHKERRA(ierr);
  }

#if 0
  // debug:= Para probar a ver si pierde memoria
  for (int tstep=1; tstep<=nstep; tstep++) {
    ierr = assemble(mesh,res,A_mom,da,xp,xpseq,x,xseq,
		    dofmap,COMP_VEC,"comp_res_mom"); CHKERRA(ierr);
  }
  PetscFinalize();
  exit(0);
#endif

  for (int tstep=1; tstep<=nstep; tstep++) {
    PetscPrintf(PETSC_COMM_WORLD,
		" --------------------------------------\n"
		"Time step: %d\n"
		" --------------------------------------\n",
		tstep);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    // FIRST (PREDICTOR) STEP

    PetscPrintf(PETSC_COMM_WORLD,"--------- Entering predictor step\n");
    // Inicializacion del paso

    ierr = VecCopy(x,dx_step);
    ierr = VecCopy(x,xold);

    for (int inwt=0; inwt<nnwt; inwt++) {
      scal=0;
      ierr = VecSet(&scal,res); CHKERRA(ierr);

      ierr = assemble(mesh,res,A_mom,da,x,xseq,xold,xoldseq,
		      dofmap,COMP_VEC,"comp_res_mom"); CHKERRA(ierr);
 
      ierr = zeroe_mat(A_mom,a_mom_assem); CHKERRA(ierr);
      
      ierr = assemble(mesh,res,A_mom,da,x,xseq,xold,xoldseq,
		      dofmap,COMP_MAT,"comp_mat_mom"); CHKERRA(ierr);
      // OJO CON LOS ARGUMENTOS X,XP AND FRIENDS. TODO ESTO HAY QUE REVISARLO
      // ierr = assemble(mesh,res,A_mom,da,xp,xpseq,x,xseq,
      // dofmap,COMP_FDJ,"comp_res_mom"); CHKERRA(ierr);

#if 0   // para ver si el jacob. analitico es igual al numerico
      ierr = assemble(mesh,res,A_mom,da,xp,xpseq,x,xseq,
		      dofmap,COMP_FDJ,"comp_res_mom"); CHKERRA(ierr);
    
      ierr = MatView(A_mom,PETSC_VIEWER_STDOUT_WORLD); CHKERRA(ierr);
#endif
    
      ierr = SLESSolve(sles_mom,res,dx,&its); CHKERRA(ierr); 
      double normres;
      ierr  = VecNorm(res,NORM_2,&normres); CHKERRA(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Newton subiter %d, norm. res-> %e\n",
		  inwt,normres); CHKERRA(ierr); 
  
      // update del subpaso
      scal= 1.0;
      ierr = VecAXPY(&scal,dx,x);
      //print_vector("outv_mom.sal",x,xseq,dofmap);
    }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    // SECOND STEP POISSON

    PetscPrintf(PETSC_COMM_WORLD,"---------- Entering poisson step\n");
    // Inicializacion del paso
    ierr = VecCopy(x,xp);

    scal=0;
    ierr = VecSet(&scal,res); CHKERRA(ierr);
    ierr = assemble(mesh,res,A_poi,da,xp,xpseq,x,xseq,
		    dofmap,COMP_VEC,"comp_res_poi"); CHKERRA(ierr);

    if (tstep==1) {
      // ierr = assemble(mesh,res,A_poi,da,xp,xpseq,x,xseq,
      // dofmap,COMP_FDJ,"comp_res_poi"); CHKERRA(ierr);
      // ierr = MatView(A_poi,matlab);

      ierr = zeroe_mat(A_poi,a_poi_assem);  CHKERRA(ierr);
      ierr = assemble(mesh,res,A_poi,da,xp,xpseq,x,xseq,
		    dofmap,COMP_MAT,"comp_mat_poi"); CHKERRA(ierr);
     }
    ierr = SLESSolve(sles_poi,res,dx,&its); CHKERRA(ierr); 

    // update del subpaso
    scal= 1.0;
    ierr = VecAXPY(&scal,dx,x);
    //    print_vector("outv_poi.sal",x,xseq,dofmap);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    // THIRD STEP PROJECTION

    PetscPrintf(PETSC_COMM_WORLD,"--------- Entering projection step\n");
    // Inicializacion del paso
    ierr = VecCopy(x,xp);

    scal=0;
    ierr = VecSet(&scal,res); CHKERRA(ierr);
    ierr = assemble(mesh,res,A_prj,da,xp,xpseq,x,xseq,
		    dofmap,COMP_VEC,"comp_res_prj"); CHKERRA(ierr);

    if (tstep==1) {
//        ierr = zeroe_mat(A_prj,a_prj_assem);  CHKERRA(ierr);
//        ierr = assemble(mesh,res,A_prj,da,xp,xpseq,x,xseq,
//  		      dofmap,COMP_FDJ,"comp_res_prj"); CHKERRA(ierr);
//        ierr = MatView(A_prj,matlab);

      ierr = zeroe_mat(A_prj,a_prj_assem);  CHKERRA(ierr);
      ierr = assemble(mesh,res,A_prj,da,xp,xpseq,x,xseq,
		      dofmap,COMP_MAT,"comp_mat_prj"); CHKERRA(ierr);
//        ierr = MatView(A_prj,matlab);

    }
    
    ierr = SLESSolve(sles_prj,res,dx,&its); CHKERRA(ierr); 

    // update del subpaso
    scal= 1.0;
    ierr = VecAXPY(&scal,dx,x);

    if (myrank==2) 	malloc_stats();
      
    // diferencia del error
    scal = -1.0;
    ierr = VecAXPY(&scal,x,dx_step);
    ierr  = VecNorm(dx_step,NORM_2,&norm); CHKERRA(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"============= ||delta u|| = %14.12e\n",norm);

    if (tstep % nsave == 0){
      print_vector("outvector.sal",x,xseq,dofmap);
    }
  }

  print_vector("outvector.sal",x,xseq,dofmap);

  ierr = SLESDestroy(sles_mom); CHKERRA(ierr);
  ierr = SLESDestroy(sles_poi); CHKERRA(ierr);
  ierr = SLESDestroy(sles_prj); CHKERRA(ierr);

  ierr = VecDestroy(x); CHKERRA(ierr); 
  ierr = VecDestroy(xseq); CHKERRA(ierr); 
  ierr = VecDestroy(xp); CHKERRA(ierr); 
  ierr = VecDestroy(xpseq); CHKERRA(ierr); 
  ierr = VecDestroy(dx); CHKERRA(ierr); 
  ierr = VecDestroy(res); CHKERRA(ierr); 

  ierr = MatDestroy(A_mom); CHKERRA(ierr); 
  ierr = MatDestroy(A_poi); CHKERRA(ierr); 
  ierr = MatDestroy(A_prj); CHKERRA(ierr); 

#endif

#ifdef NSI_TET

  Vec     x, dx, xseq, xold, xoldseq,
    xp, dx_step, xpseq,res;               /* approx solution, RHS, residual*/
  Mat     A_tet;                          /* linear system matrix */
  SLES    sles_tet;                       /* linear solver context */
  PC      pc_tet;                         /* preconditioner context */
  KSP     ksp_tet;                        /* Krylov subspace method context */
  double  norm, *sol, scal; /* norm of solution error */
  int     ierr, i, n = 10, col[3], its, flg, size, node,
    jdof, k, kk, nfixa,
    kdof, ldof, lloc, ndim, nel, nen, neq, nu,
    myrank;
  // nu:= dimension of the state vector per node
  PetscScalar  neg_one = -1.0, one = 1.0, value[3];
  PetscScalar *px;
  int jaco_num=0;
  char fcase[FLEN+1];
  Darray *da;
  Dofmap *dofmap;
  dofmap = new Dofmap;
  Mesh *mesh;

  PetscInitialize(&argc,&args,(char *)0,help);

  // Get MPI info
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);

  //  if (size != 1) SETERRA(1,0,"This is a uniprocessor example only!");
  ierr = PetscOptionsGetString(PETSC_NULL,"-case",fcase,FLEN,&flg); CHKERRA(ierr);

  // Read data

//    FILE *fid;
//    char line[LINESIZE];
//    fid=fopen(fcase,"r");
//    if (! fid) {
//      PetscPrintf(PETSC_COMM_WORLD,"Couldn't open file %s!!\n",fcase);
//      exit(1);
//    }

  // Read the mesh

  read_mesh(mesh,fcase,dofmap,neq,size,myrank,x,xseq);

  GETOPTDEF(double,tol_tet,1e-4);
  GETOPTDEF(double,tol_newton,1e-8);

  GETOPTDEF(double,atol,1e-6);
  GETOPTDEF(double,rtol,1e-3);
  GETOPTDEF(double,dtol,1e+3);
  GETOPTDEF(int,maxits,150);

  GETOPTDEF(int,nsave,10);
  GETOPTDEF(int,nstep,10000);
  GETOPTDEF(int,nnwt,1);
  GETOPTDEF(int,Krylov_dim,50);

  // initialize vectors
  ierr = VecDuplicate(x,&xp); CHKERRA(ierr);
  ierr = VecDuplicate(x,&xold); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx_step); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx); CHKERRA(ierr);
  ierr = VecDuplicate(x,&res); CHKERRA(ierr);
  ierr = VecDuplicate(xseq,&xpseq); CHKERRA(ierr);
  ierr = VecDuplicate(xseq,&xoldseq); CHKERRA(ierr);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // COMPUTE ACTVE PROFILE ON A LIBRETTO DARRAY'S
  // initialize state vectors
  scal=0;
  ierr = VecSet(&scal,x); CHKERRA(ierr);
  ierr = VecSet(&scal,xp); CHKERRA(ierr);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  char *ini_name;
  mesh->global_options->get_entry("initial_state",ini_name);
  if (ini_name!=NULL) {
    read_vector(ini_name,x,xseq,dofmap,myrank);
  } else {
    scal = 0.;
    ierr = VecSet(&scal,x); CHKERRA(ierr);
  }

  // initialize darray
  ierr = VecCopy(x,xp);
  ierr = assemble(mesh,res,A_tet,da,xp,xpseq,x,xseq,
		  dofmap,COMP_MAT_PROF,"comp_mat"); CHKERRA(ierr);

  // compute d_nnz and allocate PETSc matrix
  ierr = compute_prof(da,dofmap,myrank,&A_tet); 

  // ierr =  MatSetOption(A_tet,MAT_NEW_NONZERO_LOCATION_ERR);
  // CHKERRA(ierr);
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

  // Flags whether the matrices have been assembled or not
  int a_tet_assem=0;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // SLES para el esquema global
  ierr = SLESCreate(PETSC_COMM_WORLD,&sles_tet); CHKERRA(ierr);
  ierr = SLESSetOperators(sles_tet,A_tet,
			  A_tet,DIFFERENT_NONZERO_PATTERN); CHKERRA(ierr);
  ierr = SLESGetKSP(sles_tet,&ksp_tet); CHKERRA(ierr);
  ierr = SLESGetPC(sles_tet,&pc_tet); CHKERRA(ierr);

  ierr = KSPSetType(ksp_tet,KSPGMRES); CHKERRA(ierr);
  ierr = KSPGMRESSetRestart(ksp_tet,Krylov_dim); CHKERRA(ierr);
  ierr = KSPSetTolerances(ksp_tet,rtol,atol,dtol,maxits);
  ierr = PCSetType(pc_tet,PCJACOBI); CHKERRA(ierr);
  // ierr = PCSetType(pc_tet,PCLU); CHKERRA(ierr);
  ierr = KSPSetTolerances(ksp_tet,tol_tet,PETSC_DEFAULT,PETSC_DEFAULT,
         PETSC_DEFAULT); CHKERRA(ierr);
  ierr = KSPSetMonitor(ksp_tet,MyKSPMonitor,PETSC_NULL);


#if 0
  // debug:= Para probar a ver si pierde memoria
  for (int tstep=1; tstep<=nstep; tstep++) {
    ierr = assemble(mesh,res,A_tet,da,xp,xpseq,x,xseq,
		    dofmap,COMP_VEC,"comp_res"); CHKERRA(ierr);
  }
  PetscFinalize();
  exit(0);
#endif

#if 0   // debug:=
  PetscViewer matlab;
  ierr = ViewerASCIIOpen(PETSC_COMM_WORLD,
			 "matns.m",&matlab); CHKERRA(ierr);
  ierr = ViewerSetFormat(matlab,
			 PETSC_VIEWER_ASCII_MATLAB,"A"); CHKERRA(ierr);
  ierr = VecView(x,matlab); CHKERRA(ierr);      
#endif

  for (int tstep=1; tstep<=nstep; tstep++) {
    PetscPrintf(PETSC_COMM_WORLD,
		" --------------------------------------\n"
		"Time step: %d\n",
		tstep);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    // TET ALGORITHM

    // Inicializacion del paso

    ierr = VecCopy(x,dx_step);
    ierr = VecCopy(x,xold);

    for (int inwt=0; inwt<nnwt; inwt++) {
      scal=0;
      ierr = VecSet(&scal,res); CHKERRA(ierr);

      ierr = VecCopy(x,xp);
      ierr = assemble(mesh,res,A_tet,da,x,xseq,xold,xoldseq,
		      dofmap,COMP_VEC,"comp_res"); CHKERRA(ierr);

      //ierr = zeroe_mat(A_tet,a_tet_assem); CHKERRA(ierr);
      ierr = MatZeroEntries(A_tet); CHKERRA(ierr);
      ierr = assemble(mesh,res,A_tet,da,x,xseq,xold,xoldseq,
		      dofmap,COMP_MAT,"comp_mat"); CHKERRA(ierr);

      // if (myrank==2) malloc_stats();
      
      //      show_mallinfo();
      if (tstep==1) 
	ierr =  MatSetOption(A_tet,MAT_NO_NEW_NONZERO_LOCATIONS);

#if 0
      scal=-1;
      ierr = MatScale(&scal,A_tet);
      ierr = assemble(mesh,res,A_tet,da,xp,xpseq,x,xseq,
		      dofmap,COMP_FDJ,"comp_res"); CHKERRA(ierr);
      ierr = MatView(A_tet,matlab); CHKERRA(ierr);
            PetscFinalize();
      exit(0);
#endif
      
      ierr = SLESSolve(sles_tet,res,dx,&its); CHKERRA(ierr); 
      double normres;
      ierr  = VecNorm(res,NORM_2,&normres); CHKERRA(ierr);
  
      // update del subpaso
      scal= 1.0;
      ierr = VecAXPY(&scal,dx,x);
    }

    // diferencia del error
    scal = -1.0;
    ierr = VecAXPY(&scal,x,dx_step);
    ierr  = VecNorm(dx_step,NORM_2,&norm); CHKERRA(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"============= ||delta u|| = %14.12e\n",norm);
    if (norm<tol_newton) break;

    if (tstep % nsave == 0){
#if 0
      PetscViewer matlab;
      ierr = ViewerASCIIOpen(PETSC_COMM_WORLD,
			     "mat.output",&matlab); CHKERRA(ierr);
      ierr = ViewerSetFormat(matlab,
      			     PETSC_VIEWER_ASCII_MATLAB,"pepe"); CHKERRA(ierr);
      ierr = VecView(x,matlab); CHKERRA(ierr);      
      ViewerDestroy(matlab);
#endif 
      //      ierr = VecView(res,PETSC_VIEWER_STDOUT_WORLD); CHKERRA(ierr);      
      //PetscFinalize();
      //exit(0);
      print_vector("outvector.sal",x,xseq,dofmap);
    }
  }

  print_vector("outvector.sal",x,xseq,dofmap);

  ierr = SLESDestroy(sles_tet); CHKERRA(ierr);

  ierr = VecDestroy(x); CHKERRA(ierr); 
  ierr = VecDestroy(xseq); CHKERRA(ierr); 
  ierr = VecDestroy(xp); CHKERRA(ierr); 
  ierr = VecDestroy(xpseq); CHKERRA(ierr); 
  ierr = VecDestroy(dx); CHKERRA(ierr); 
  ierr = VecDestroy(res); CHKERRA(ierr); 

  ierr = MatDestroy(A_tet); CHKERRA(ierr); 

#endif

  PetscFinalize();
  exit(0);
}
