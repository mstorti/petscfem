//__INSERT_LICENSE__
// $Id: dvecpar.cpp,v 1.5 2006/04/08 20:37:15 mstorti Exp $

#include <src/dvecpar.h>

#define DVECTOR_MPI_TYPE_DEF(stype,TYPE)			\
MPI_Datatype							\
dvector_mpi_type<stype>::type() { return MPI_##TYPE; } ;	\
								\
// Instantiates for this type					\
template							\
void								\
dvector_clone_parallel(dvector<type> &w,int root);		\
								\
template							\
void								\
dvector_read_parallel(const char *file,				\
		      dvector<type> &w,int root=0);


DVECTOR_MPI_TYPE_DEF(int,INT) 
DVECTOR_MPI_TYPE_DEF(double,DOUBLE) 
DVECTOR_MPI_TYPE_DEF(char,CHAR) 
DVECTOR_MPI_TYPE_DEF(float,FLOAT) 
