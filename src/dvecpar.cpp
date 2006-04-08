//__INSERT_LICENSE__
// $Id: dvecpar.cpp,v 1.6 2006/04/08 22:29:32 mstorti Exp $

#include <src/dvecpar.h>
#include <src/dvecpar2.h>

#define DVECTOR_MPI_TYPE_DEF(stype,TYPE)			\
MPI_Datatype							\
dvector_mpi_type<stype>::type() { return MPI_##TYPE; } ;	\

#if 0
// Instantiates for this type					\
template							\
void								\
dvector_clone_parallel(dvector<type> &w,int root);		\
								\
template							\
void								\
dvector_read_parallel(const char *file,				\
		      dvector<type> &w,int root=0);
#endif

#if 0
DVECTOR_MPI_TYPE_DEF(int,INT) 
DVECTOR_MPI_TYPE_DEF(double,DOUBLE) 
DVECTOR_MPI_TYPE_DEF(char,CHAR) 
DVECTOR_MPI_TYPE_DEF(float,FLOAT) 
#endif

MPI_Datatype						       
dvector_mpi_type<int>::type() { return MPI_INT; } ;	
