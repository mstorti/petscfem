//__INSERT_LICENSE__
// $Id: dvecpar.cpp,v 1.1 2006/04/06 16:33:12 mstorti Exp $

#include <src/dvecpar.h>

MPI_Datatype dvector_mpi_type(int) { return MPI_INT; }
MPI_Datatype dvector_mpi_type(double) { return MPI_DOUBLE; }
MPI_Datatype dvector_mpi_type(char) { return MPI_CHAR; }
