/*      "$Id: try.cpp,v 1.1.2.1 2004/07/07 14:31:30 mstorti Exp $"; */

static char help[] = "Demonstrates using ISLocalToGlobalMappings.\n\n";

/*T
  Concepts: local to global mappings
  Concepts: global to local mappings
  
  Description:  Creates an index set based on blocks of integers. Views that index set
  and then destroys it.
  T*/

#include "petscis.h"

int main(int argc,char **argv)
{
  ISLocalToGlobalMapping mapping;
  
  int myrank,size;
  PetscInitialize(&argc,&argv,NULL,NULL);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
  //  int                    i,n = 4,ierr,indices[] = {0,3,9,12},m = 2,input[] = {0,2};
  int output[2],inglobals[13],outlocals[13];
  int i,n = 4,ierr,indices[] = {myrank,myrank+3,myrank+9,myrank+12},m = 2,input[] = {0,2};
  /*
    Create a local to global mapping. Each processor independently
    creates a mapping  
  */
  ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD,n,indices,&mapping);
  
  /*
    Map a set of local indices to their global values 
  */
  ISLocalToGlobalMappingApply(mapping,m,input,output);
  PetscIntView(n,indices,PETSC_VIEWER_STDOUT_SELF);
  PetscIntView(m,output,PETSC_VIEWER_STDOUT_SELF);
  
  /*
    Map some global indices to local, retaining the ones without a local index by -1
  */
  for (i=0; i<13; i++) {
    inglobals[i] = i;
  }
  ISGlobalToLocalMappingApply(mapping,IS_GTOLM_MASK,13,inglobals,PETSC_NULL,outlocals);
  
  PetscIntView(13,outlocals,PETSC_VIEWER_STDOUT_SELF);
  
  /*
    Map some global indices to local, dropping the ones without a local index.
  */
  ISGlobalToLocalMappingApply(mapping,IS_GTOLM_DROP,13,inglobals,&m,outlocals);
  
  PetscIntView(m,outlocals,PETSC_VIEWER_STDOUT_SELF);
  
  /*
    Free the space used by the local to global mapping
  */
  ISLocalToGlobalMappingDestroy(mapping);
  
  
  PetscFinalize();
  return 0;
}


