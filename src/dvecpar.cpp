//__INSERT_LICENSE__
// $Id: dvecpar.cpp,v 1.2 2006/04/06 16:37:44 mstorti Exp $

#include <src/dvecpar.h>

#define DVECTOR_MPI_TYPE_DEF(type,TYPE) 
MPI_Datatype dvector_mpi_type(type) { return MPI_##TYPE; }

DVECTOR_MPI_TYPE_DEF(int,INT) 
DVECTOR_MPI_TYPE_DEF(double,DOUBLE) 
DVECTOR_MPI_TYPE_DEF(char,CHAR) 
DVECTOR_MPI_TYPE_DEF(float,FLOAT) 

#if 0
MPI_Datatype dvector_mpi_type(int) { return MPI_INT; }
MPI_Datatype dvector_mpi_type(double) { return MPI_DOUBLE; }
MPI_Datatype dvector_mpi_type(char) { return MPI_CHAR; }
MPI_Datatype dvector_mpi_type(float) { return MPI_FLOAT; }
#endif
