/*$Id: main.c,v 1.1 2004/07/07 17:01:44 rodrigop Exp $*/
static char help[] = "Solves 2d-laplacian on quadrilateral grid.\n\
   Options:\n\
    -show_solution pipe solution to matlab (visualized with bscript.m).\n\
    -show_griddata print the local index sets and local to global mappings \n\
    -show_matrix visualize the sparsity structure of the stiffness matrix.\n\
    -show_grid visualize the global and local grids with numbering.\n";

/*
    The file appctx.h includes all the data structures used by this code
*/
#include "appctx.h"


EXTERN_C_BEGIN
extern int PCCreate_NN (PC);
extern int MatPartitioningCreate_Square (MatPartitioning);
EXTERN_C_END

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  int            ierr,its,rank;
  AppCtx         *appctx;     /* contains all the data used by this PDE solver */

  /* ---------------------------------------------------------------------
     Initialize PETSc
     --------------------------------------------------------------------- */

  PetscFunctionBegin;

  PetscInitialize(&argc,&argv,PETSC_NULL,help);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  PCRegisterDynamic("nn",PETSC_NULL,"PCCreate_NN",PCCreate_NN);
  MatPartitioningRegisterDynamic("square",PETSC_NULL,"MatPartitioningCreate_Square",MatPartitioningCreate_Square);
                                                      
  /*  Load the grid database -- in appload.c              */
  ierr = AppCtxCreate(PETSC_COMM_WORLD,&appctx);CHKERRQ(ierr);

  /*   Setup the graphics routines to view the grid -- in appview.c  */
  ierr = AppCtxGraphics(appctx);CHKERRQ(ierr);

 
  /*   Setup the linear system and solve it -- in appalgebra.c */
  ierr = AppCtxSolve(appctx,&its);CHKERRQ(ierr);
  if (appctx->view.save_matrix) {  
    PetscViewer ppview;
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"advdiff_lin_matrix.m",&ppview);CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(ppview,PETSC_VIEWER_ASCII_MATLAB);
    ierr = MatView(appctx->algebra.A,ppview);CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(ppview,PETSC_VIEWER_ASCII_MATLAB);    
    ierr = VecView(appctx->algebra.b,ppview);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(ppview);CHKERRQ(ierr);
  }

  /*   Send solution to  matlab viewer -- in appview.c */
  if (appctx->view.show_solution) {
    ierr = AppCtxViewMatlab(appctx);CHKERRQ(ierr);  
  }

  /*  Destroy all datastructures  -- in appload.c */
  ierr = AppCtxDestroy(appctx);CHKERRQ(ierr);

  /* Close down PETSc and stop the program */
  ierr = PetscFinalize();CHKERRQ(ierr);

  // if (!rank) { printf("\n\nNumber of Iterations: %d\n\n",its); }

  PetscFunctionReturn(0);
}










