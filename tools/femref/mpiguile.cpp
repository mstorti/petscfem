#include <cassert>
#include <cstdio>
#include <cmath>
#include <libguile.h>
#if 0
#include "vector.h"
#include <petsc.h>
#include "./petscscm.h"
#include "./dvector.h"
#endif
#include <mpi.h>

#define N 5

typedef SCM(*scm_fun)();

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ "mpi-send"
static SCM
mpi_send(SCM s_val,SCM s_dest) {
  SCM_ASSERT(scm_number_p(s_val),s_val,SCM_ARG1,__FUN__);
  double val = scm_num2double(s_val,0,__FUN__);
  double tval = round(val);

  SCM_ASSERT(SCM_INUMP(s_dest),
	     s_dest,SCM_ARG2,__FUN__);
  int dest = SCM_INUM(s_dest);

  printf("mpi_send: sending %lg, error %lg\n",val,val-tval);
  double v[N];
  for (int j=0; j<N; j++) v[j] = val;
  printf("mpi_send: sending %lg, error %lg\n",v[0],v[0]-tval);
  MPI_Send(v,N,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
  return SCM_UNSPECIFIED;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ "mpi-recv"
static SCM
mpi_recv(SCM s_source) {
  SCM_ASSERT(SCM_INUMP(s_source),
	     s_source,SCM_ARG1,__FUN__);
  int source = SCM_INUM(s_source);
  double v[N];
#if 0
  for (int j=0; j<N; j++)
    v[j] = 0.1234567890123456789;
  for (int j=0; j<N; j++)
    printf("%g ",v[j]);
  printf(" xxx \n");
#endif
  MPI_Status status;
  MPI_Recv(v,N,MPI_DOUBLE,source,0,
	   MPI_COMM_WORLD,&status);
  double val = v[0];
  double tval = round(val);
  printf("mpi_recv error: received ");
  for (int j=0; j<N; j++)
    printf("%g ",v[j]-tval);
  printf("\n");
#if 0
  for (int j=1; j<N; j++)
    assert(v[j]==val);
#endif
  return scm_make_real(val);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ "mpi-rank"
static SCM
mpi_rank() {
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  return SCM_MAKINUM(myrank);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ "mpi-size"
static SCM
mpi_size() {
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  return SCM_MAKINUM(size);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ "mpi-finalize"
static SCM
mpi_finalize() {
  MPI_Finalize();
  return SCM_UNSPECIFIED;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
extern "C" void
init_mpi (void) {
  scm_c_define_gsubr("mpi-send",2,0,0,scm_fun(mpi_send));
  scm_c_define_gsubr("mpi-recv",1,0,0,scm_fun(mpi_recv));
  scm_c_define_gsubr("mpi-rank",0,0,0,scm_fun(mpi_rank));
  scm_c_define_gsubr("mpi-size",0,0,0,scm_fun(mpi_size));
  scm_c_define_gsubr("mpi-finalize",0,0,0,scm_fun(mpi_finalize));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
static void
inner_main (void *closure, int argc, char **argv) {
  init_mpi();
  scm_shell(argc, argv);
  MPI_Finalize();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int main (int argc, char **argv) {
  MPI_Init(&argc,&argv);
  scm_boot_guile (argc, argv, inner_main, 0);
  return 0; // never reached
}

