//__INSERT_LICENSE__
// $Id: dvecpar.cpp,v 1.3 2006/04/06 16:38:14 mstorti Exp $

#include <src/dvecpar.h>

#define DVECTOR_MPI_TYPE_DEF(type,TYPE) 
MPI_Datatype dvector_mpi_type(type) { return MPI_##TYPE; }

DVECTOR_MPI_TYPE_DEF(int,INT) 
DVECTOR_MPI_TYPE_DEF(double,DOUBLE) 
DVECTOR_MPI_TYPE_DEF(char,CHAR) 
DVECTOR_MPI_TYPE_DEF(float,FLOAT) 
