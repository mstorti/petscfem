/*$Id: appalgebra.c,v 1.1.2.1 2004/03/01 21:54:57 mstorti Exp $*/
#include "appctx.h"
#include "math.h"
#include "chrono.h"
//#include "petscmat.h"
extern int rank;

/*
         Sample right hand side and boundary conditions
*/
#undef __FUNCT__
#define __FUNCT__ "pde_rhs"
int pde_rhs(void *dummy,int n,double *xx,double *f)
{
  double pi = PETSC_PI, x = xx[0], y = xx[1];
  PetscFunctionBegin;
  // *f = 8*pi*pi*sin(2*pi*x)*sin(2*pi*y)-20*pi*cos(2*pi*x)*sin(2*pi*y);
  *f = 1.0;
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "pde_bc"
int pde_bc(void *dummy,int n,double *xx,double *f)
{
  double pi = 3.1415927, x = xx[0], y = xx[1];
  PetscFunctionBegin;
  //  *f = sin(2*pi*x)*sin(2*pi*y);
  *f = 0.0;
  PetscFunctionReturn(0);
}

/*
  Sets up the linear system associated with the PDE and solves it
*/
#undef __FUNCT__
#define __FUNCT__ "AppCxtSolve"
int AppCtxSolve(AppCtx* appctx, int *its)
{
  AppAlgebra  *algebra = &appctx->algebra;
  MPI_Comm    comm = appctx->comm;
  SLES        sles;
  int         ierr,myrank;
  Chrono chrono;
  PetscFunctionBegin;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);CHKERRQ(ierr);
  /*  Set the functions to use for the right hand side and Dirichlet boundary */
  ierr = PFCreate(comm,2,1,&appctx->element.rhs);CHKERRQ(ierr);
  ierr = PFSetOptionsPrefix(appctx->element.rhs,"rhs_");CHKERRQ(ierr);
  ierr = PFSetType(appctx->element.rhs,PFQUICK,(void*)pde_rhs);CHKERRQ(ierr);
  ierr = PFSetFromOptions(appctx->element.rhs);CHKERRQ(ierr);

  ierr = PFCreate(comm,2,1,&appctx->bc);CHKERRQ(ierr);
  ierr = PFSetOptionsPrefix(appctx->bc,"bc_");CHKERRQ(ierr);
  ierr = PFSetType(appctx->bc,PFQUICK,(void*)pde_bc);CHKERRQ(ierr);
  ierr = PFSetFromOptions(appctx->bc);CHKERRQ(ierr);

  /*     A) Set the quadrature values for the reference element  */
  ierr = SetReferenceElement(appctx);CHKERRQ(ierr);

  /*     1) Create vector to contain load and various work vectors  */
  ierr = AppCtxCreateRhs(appctx);CHKERRQ(ierr);

  /*     2)  Create the sparse matrix,with correct nonzero pattern  */
  ierr = AppCtxCreateMatrix(appctx);CHKERRQ(ierr);

  /*     3)  Set the right hand side values into the load vector   */
  ierr = AppCtxSetRhs(appctx);CHKERRQ(ierr);

  /*     4)  Set the matrix entries   */
  chrono.start();
  double t_ass=chrono.elapsed();
  ierr = AppCtxSetMatrix(appctx);CHKERRQ(ierr);
  double t_ass_f=chrono.elapsed()-t_ass;
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			  "Time spent in Setting Matrix[proc %d]:  %f,\n",myrank,t_ass_f);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  /* view sparsity structure of the matrix */
  if (appctx->view.show_matrix) {  
    ierr = PetscPrintf(PETSC_COMM_WORLD,"The stiffness matrix, before bc applied\n");CHKERRQ(ierr);
    ierr = MatView(appctx->algebra.A,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);
  }

  /*     6) Set the matrix boundary conditions */
  ierr = SetMatrixBoundaryConditions(appctx);CHKERRQ(ierr);

  /* view sparsity structure of the matrix */
  if(appctx->view.show_matrix) {  
    ierr = PetscPrintf(PETSC_COMM_WORLD,"The stiffness matrix, after bc applied\n");CHKERRQ(ierr);
    ierr = MatView(appctx->algebra.A,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);
  }

  PreLoadBegin(PETSC_TRUE,"Solver setup");  

    /*     5) Set the rhs boundary conditions - this also creates initial guess that satisfies boundary conditions */
    ierr = SetBoundaryConditions(appctx);CHKERRQ(ierr);

    ierr = SLESCreate(comm,&sles);CHKERRQ(ierr);
    ierr = SLESSetOperators(sles,algebra->A,algebra->A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    KSP ksp;
    ierr = SLESGetKSP(sles,&ksp);CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
    double rtol=1.0e-9,atol=1.0e-14,dtol=1.0e5;
    int maxits=200;
    
    ierr = KSPSetTolerances(ksp,rtol,atol,dtol,maxits); CHKERRQ(ierr); 
    ierr = KSPSetMonitor(ksp,KSPDefaultMonitor,NULL,NULL); CHKERRQ(ierr); 
    ierr = KSPSetType(ksp,KSPGMRES);
    ierr = KSPGMRESSetRestart(ksp,maxits);
    //  ierr = KSPGMRESSetOrthogonalization(ksp,KSPGMRESIROrthogonalization);
    ierr = KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization);

    ierr = SLESSetFromOptions(sles);CHKERRQ(ierr);
    ierr = SLESSetUp(sles,appctx->algebra.b,appctx->algebra.b);CHKERRQ(ierr);
    ierr = SLESSetFromOptions(sles);CHKERRQ(ierr);
    ierr = SLESSetUp(sles,appctx->algebra.b,appctx->algebra.b);CHKERRQ(ierr);
    
    PreLoadStage("Solve");  
    double t_ini=chrono.elapsed();
    ierr = SLESSolve(sles,algebra->b,algebra->x,its);CHKERRQ(ierr);
    double t_end=chrono.elapsed()-t_ini;
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			    "Time spent in Solving System[proc %d]:  %f,\n",myrank,t_end);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    {
      PetscTruth flg;
      ierr = PetscOptionsHasName(PETSC_NULL,"-save_global_preconditioner",&flg);CHKERRQ(ierr);
      if (flg) {
	PC          pc;
	KSP         ksp;
	Mat         mat,mat2;
	PetscViewer viewer;
	ierr = SLESGetPC(sles,&pc);CHKERRQ(ierr);
	ierr = SLESGetKSP(sles,&ksp);CHKERRQ(ierr);
	ierr = PCComputeExplicitOperator(pc,&mat);CHKERRQ(ierr);
	ierr = KSPComputeExplicitOperator(ksp,&mat2);CHKERRQ(ierr);
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"pc.m",&viewer);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
	ierr = MatView(mat,viewer);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
	ierr = MatView(mat2,viewer);CHKERRQ(ierr);
	ierr = MatDestroy(mat);CHKERRQ(ierr);
	ierr = MatDestroy(mat2);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
      }
    }

    /*      Free the solver data structures */
    ierr = SLESDestroy(sles);CHKERRQ(ierr);
  PreLoadEnd();

  ierr = PFDestroy(appctx->bc);CHKERRQ(ierr);
  ierr = PFDestroy(appctx->element.rhs);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*----------------------------------------------------------------
       1  -  Generates the "global" parallel vector to contain the 
	     right hand side and solution.
*/
#undef __FUNCT__
#define __FUNCT__ "AppCxtCreateRhs"
int AppCtxCreateRhs(AppCtx *appctx)
{
  AppGrid     *grid = &appctx->grid;
  AppAlgebra  *algebra = &appctx->algebra;
  MPI_Comm    comm = appctx->comm;
  int         ierr;

  PetscFunctionBegin;
  /*  Create vector to contain load,  local size should be number of  vertices  on this proc.  */
  ierr = VecCreateMPI(comm,grid->vertex_local_n,PETSC_DETERMINE,&algebra->b);CHKERRQ(ierr);

  /* This allows one to set entries into the vector using the LOCAL numbering: via VecSetValuesLocal() */
  ierr = VecSetLocalToGlobalMapping(algebra->b,grid->ltog);CHKERRQ(ierr);

  /* Generate the vector to contain the solution */
  ierr = VecDuplicate(algebra->b,&algebra->x);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*---------------------------------------------------------------
      2  - Generates the "global" parallel matrix
*/
#undef __FUNCT__
#define __FUNCT__ "AppCxtCreateMatrix"
int AppCtxCreateMatrix(AppCtx* appctx)
{

  AppAlgebra  *algebra = &appctx->algebra;
  AppGrid     *grid    = &appctx->grid;
  MPI_Comm    comm = appctx->comm;
  int         ierr;
 
  PetscFunctionBegin;
  ierr = MatCreate(comm,grid->vertex_local_n,grid->vertex_local_n,PETSC_DETERMINE,PETSC_DETERMINE,&algebra->A);CHKERRQ(ierr);
  ierr = MatSetFromOptions(algebra->A);CHKERRQ(ierr);

  /* Allows one to set values into the matrix using the LOCAL numbering, via MatSetValuesLocal() */
  ierr = MatSetLocalToGlobalMapping(algebra->A,grid->ltog);  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*---------------------------------------------------------------------------
     3 - Computes the entries in the right hand side and sets them into the parallel vector
         Uses B and C
*/
#undef __FUNCT__
#define __FUNCT__ "AppCxtSetRhs"
int AppCtxSetRhs(AppCtx* appctx)
{
  /********* Context informatrion ***********/
  AppGrid    *grid = &appctx->grid;
  AppAlgebra *algebra = &appctx->algebra;
  AppElement *phi = &appctx->element;

  /****** Local Variables ***********/
  int        ierr,i;
  int        *vertex_ptr;
  int        bn = 4; /* number of basis functions */
  int        vertexn = 4; /* number of degrees of freedom  */

  PetscFunctionBegin;
  /* loop over local cells */
  for(i=0;i<grid->cell_n;i++){

    /* coords_ptr points to the coordinates of the current cell */
    phi->coords = grid->cell_coords + 2*bn*i;  /*number of cell coords */

    /* compute the values of basis functions on this element */
    ierr = SetLocalElement(phi);CHKERRQ(ierr); 

    /* compute the  element load (integral of f with the 4 basis elements)  */
    /* values get put into phi->rhsresult  */
    ierr = ComputeRHSElement(phi);CHKERRQ(ierr);

    /*********  Set Values *************/
    /* vertex_ptr points to place in the vector to set the values */
    vertex_ptr = grid->cell_vertex + vertexn*i; 

    ierr = VecSetValuesLocal(algebra->b,bn,vertex_ptr,phi->rhsresult,ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(algebra->b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(algebra->b);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}  

/*------------------------------------------------------------------
      4 - Computes the element stiffness matrices and stick into 
   global stiffness matrix. Uses B and D.
*/
#undef __FUNCT__
#define __FUNCT__ "AppCxtSetMatrix"
int AppCtxSetMatrix(AppCtx* appctx)
{
  /********* Contex information ***********/
  AppAlgebra *algebra = &appctx->algebra;
  AppGrid    *grid    = &appctx->grid;
  AppElement *phi = &appctx->element; 

  /****** Local Variables ***********/
  int        i,ierr;
  int        *vertex_ptr;
  int        bn = 4; /* number of basis functions */
  int        vertexn = 4; /* number of degrees of freedom */
  Mat local;//for preallocation
  PetscFunctionBegin;

  //  ierr = MatISGetLocalMat(algebra->A,&local);CHKERRQ(ierr);
  //  ierr = MatSeqAIJSetPreallocation(local,9,PETSC_NULL);CHKERRQ(ierr);

  /* loop over local cells */
  for(i=0;i<grid->cell_n;i++){

    /* coords_ptr points to the coordinates of the current cell */
    phi->coords = grid->cell_coords + 2*bn*i;/*number of cell coords */

    /* compute the values of basis functions on this element */
    ierr = SetLocalElement(phi);CHKERRQ(ierr);
   
    /*    Compute the element stiffness  */  
    /* result is returned in phi->stiffnessresult */
    ierr = ComputeStiffnessElement(phi);CHKERRQ(ierr);

    /*********  Set Values *************/
    /* vertex_ptr points to place in the matrix to set the values */
    vertex_ptr = grid->cell_vertex + vertexn*i;

    ierr = MatSetValuesLocal(algebra->A,vertexn,vertex_ptr,vertexn,vertex_ptr,(double*)phi->stiffnessresult,ADD_VALUES);CHKERRQ(ierr);
    /* ierr = MatSetValues(algebra->localA,vertexn,vertex_ptr,vertexn,vertex_ptr,(double*)phi->stiffnessresult,ADD_VALUES);CHKERRQ(ierr);*/
  }
  ierr = MatAssemblyBegin(algebra->A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(algebra->A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*----------------------------------------------------------------
      5   - Apply the Dirichlet boundary conditions (see 6 also).
     This places the Dirichlet function value on the right hand side
     and 6 sticks a row of the identity matrix on the left side 
     thus forcing the solution at the given points to match the 
     Dirichlet function.
*/
#undef __FUNCT__
#define __FUNCT__ "SetBoundaryConditions"
int SetBoundaryConditions(AppCtx *appctx)
{
 /********* Context informatrion ***********/
  AppAlgebra *algebra = &appctx->algebra;
  AppGrid    *grid = &appctx->grid;

  /****** Local Variables ***********/
  int          ierr,i;
  int          *vertex_ptr; 
  PetscScalar  zero = 0.0;

  PetscFunctionBegin;

  /*  -------------------------------------------------------------
         Apply Dirichlet boundary conditions
      -----------------------------------------------------------*/
  
  /* need to set the points on RHS corresponding to vertices on the boundary to
     the desired value. */

  /* get list of vertices on the boundary */
  ierr = ISGetIndices(grid->vertex_boundary,&vertex_ptr);CHKERRQ(ierr);
  for(i=0;i<grid->boundary_n;i++){
    /* evaluate boundary condition function at point */
    ierr = PFApply(appctx->bc,1,&grid->boundary_coords[2*i],&grid->boundary_values[i]);CHKERRQ(ierr);
  }

  /* set the right hand side values at those points */
  ierr = VecSetValuesLocal(algebra->b,grid->boundary_n,vertex_ptr,grid->boundary_values,INSERT_VALUES);CHKERRQ(ierr);

  /* set initial guess satisfying boundary conditions */
  ierr = VecSet(&zero,algebra->x);CHKERRQ(ierr);
  ierr = VecSetValuesLocal(algebra->x,grid->boundary_n,vertex_ptr,grid->boundary_values,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(algebra->x);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(algebra->x);CHKERRQ(ierr);

  ierr = ISRestoreIndices(grid->vertex_boundary,&vertex_ptr);CHKERRQ(ierr);
 
  ierr = VecAssemblyBegin(algebra->b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(algebra->b);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------
     6 - Set the matrix boundary conditions (see also 5). Replace the corresponding 
         rows in the matrix with the identity.
*/
#undef __FUNCT__
#define __FUNCT__ "SetMatrixBoundaryConditions"
int SetMatrixBoundaryConditions(AppCtx *appctx)
{
  /********* Context informatrion ***********/
  AppAlgebra *algebra = &appctx->algebra;
  AppGrid    *grid = &appctx->grid;

  /****** Local Variables ***********/
  double     one = 1.0;
  int        ierr;

  PetscFunctionBegin;
  ierr = MatZeroRowsLocal(algebra->A,grid->vertex_boundary,&one);CHKERRQ(ierr); 
  PetscFunctionReturn(0);
}
