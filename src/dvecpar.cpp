//__INSERT_LICENSE__
// $Id: dvecpar.cpp,v 1.9 2007/01/30 19:03:44 mstorti Exp $

#include <src/dvecpar.h>
#include <src/dvecpar2.h>

#define DVECTOR_MPI_TYPE_DEF(stype,TYPE)			\
template<> MPI_Datatype						\
dvector_mpi_type<stype>::type() { return MPI_##TYPE; }  	\
template<>							\
void								\
dvector_clone_parallel(dvector<stype> &w,int root);		\
								\
template<>							\
void								\
dvector_read_parallel(const char *file,				\
		      dvector<stype> &w,int root);

DVECTOR_MPI_TYPE_DEF(int,INT) 
DVECTOR_MPI_TYPE_DEF(double,DOUBLE) 
DVECTOR_MPI_TYPE_DEF(char,CHAR) 
DVECTOR_MPI_TYPE_DEF(float,FLOAT) 
