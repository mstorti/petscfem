/*__INSERT_LICENSE__*/
// $Id: pfmat.cpp,v 1.1.2.10 2001/12/29 23:59:35 mstorti Exp $

// Tests for the `PFMat' class
#include <src/debug.h>
#include <src/fem.h>
#include <src/utils.h>
#include <src/elemset.h>
#include <src/pfmat.h>
#include <src/pfptscmat.h>
#include <src/iisdmat.h>
#include <src/petscmat.h>
#include <src/graph.h>

extern int MY_RANK,SIZE;

class Part : public DofPartitioner {
public:
  int N,comm_size;
  int processor(int k) const { return k*comm_size/N; }
  ~Part() {}
} part;

#define CHKOPT(name)  assert(!strcmp(#name,args[++arg]))

int main(int argc,char **args) {

  PFMat *A_p;
  Vec x,b,xex;
  int ierr,iisd_subpart,nsolve,rand_rhs,mat_type;
  // debug.activate();
  int myrank,size;
  PetscInitialize(&argc,&args,NULL,NULL);

  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);

  double cond,L,Q,rand_coef=1.;
  int Nelem=10, debug_print=0, nmat;
  int arg=0;

  // usage: pfmat.bin ne <Nelem> deb <debug_print> nsl <nsolve> sbdp
  // <iisd_subpart> nmat <nmat> rrhs <rand_rhs> mtyp <mat_type>

  CHKOPT(ne);
  if (myrank==0 && argc>++arg)
    sscanf(args[arg],"%d",&Nelem);
  ierr = MPI_Bcast (&Nelem, 
		    1, MPI_INT, 0,MPI_COMM_WORLD); CHKERRA(ierr); 
  int N=Nelem-1;

  CHKOPT(deb);
  if (myrank==0 && argc>++arg)
    sscanf(args[arg],"%d",&debug_print);
  ierr = MPI_Bcast (&debug_print, 
		    1, MPI_INT, 0,MPI_COMM_WORLD); CHKERRA(ierr); 

  CHKOPT(nsl);
  nsolve=1;
  if (myrank==0 && argc>++arg)
    sscanf(args[arg],"%d",&nsolve);
  ierr = MPI_Bcast (&nsolve, 1, 
		    MPI_INT, 0,MPI_COMM_WORLD); CHKERRA(ierr); 
    
  CHKOPT(sbdp);
  iisd_subpart=1;
  if (myrank==0 && argc>++arg) 
    sscanf(args[arg],"%d",&iisd_subpart);
  ierr = MPI_Bcast (&iisd_subpart, 
		    1, MPI_INT, 0,MPI_COMM_WORLD); CHKERRA(ierr); 

  CHKOPT(nmat);
  nmat=1;
  if (myrank==0 && argc>++arg) 
    sscanf(args[arg],"%d",&nmat);
  ierr = MPI_Bcast (&nmat, 
		    1, MPI_INT, 0,MPI_COMM_WORLD); CHKERRA(ierr); 

  CHKOPT(rrhs);
  rand_rhs=0;
  if (myrank==0 && argc>++arg) 
    sscanf(args[arg],"%d",&rand_rhs);
  ierr = MPI_Bcast (&rand_rhs, 
		    1, MPI_INT, 0,MPI_COMM_WORLD); CHKERRA(ierr); 

  CHKOPT(mtyp);
  mat_type=0; // 0 -> IISDMat, 1 -> PETScMat
  if (myrank==0 && argc>++arg) 
    sscanf(args[arg],"%d",&mat_type);
  ierr = MPI_Bcast (&mat_type, 
		    1, MPI_INT, 0,MPI_COMM_WORLD); CHKERRA(ierr); 

  PetscPrintf(PETSC_COMM_WORLD,
	      "iisd_subpart %d, nsolve %d, debug_print %d, Nelem %d\n",
	      iisd_subpart, nsolve, debug_print, Nelem);

  part.comm_size = size;
  part.N = N;
  MY_RANK = myrank;
  SIZE = size;

  IISDMat AA(N,N,part,PETSC_COMM_WORLD);
  PETScMat AAA(N,N,part,PETSC_COMM_WORLD);
  if (mat_type==0) {
    A_p = &AA;
  } else if (mat_type==0) {
    A_p = &AAA;
  } else assert(0);
  PFMat &A = *A_p;

  for (int j=0; j<N; j++) {
    if (j % size != myrank) continue; // Load periodically
    A.set_profile(j,j);
    if (j+1<N) A.set_profile(j,j+1);
    if (j-1>=0) A.set_profile(j,j-1);
  }
  A.set_option("iisd_subpart",iisd_subpart);
  A.set_option("print_internal_loop_conv",debug_print);
  A.set_option("rtol",1e-8);
  A.set_option("atol",0);
  A.create();
  if (debug_print) A.view();

  int nhere=0;
  for (int k=0; k<N; k++) 
    if (part.processor(k)==myrank) nhere++;
  
  ierr = VecCreateMPI(PETSC_COMM_WORLD,
		      nhere,PETSC_DETERMINE,&b); CHKERRA(ierr); 
  ierr = VecDuplicate(b,&x); CHKERRA(ierr); 
  ierr = VecDuplicate(b,&xex); CHKERRA(ierr); 

  for (int imat=0; imat<nmat; imat++) {
    A.zero_entries();
    if (myrank==0) {
      cond = 1. + rand_coef * ::drand();
      L = 1. + rand_coef * ::drand();
    }

    ierr = MPI_Bcast (&cond, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    ierr = MPI_Bcast (&L, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    double h = L/double(Nelem);
    double coef = cond / (h*h);

    for (int j=0; j<N; j++) {
      if (j % size != myrank) continue; // Load periodically
      A.set_value(j,j,2.*coef);
      if (j+1 <  N ) A.set_value(j,j+1,-coef);
      if (j-1 >= 0 ) A.set_value(j,j-1,-coef);
    }
    A.assembly_begin(MAT_FINAL_ASSEMBLY); 
    A.assembly_end(MAT_FINAL_ASSEMBLY); 

    for (int ksolve=0; ksolve<nsolve; ksolve++) {

      if (myrank==0) {
	Q = 1. + rand_coef * ::drand();
	printf("cond %f, Q %f, L %f\n",cond,Q,L);
      }
      ierr = MPI_Bcast (&Q, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);

      if (myrank==0) {
	for (int j=0; j<N; j++) {
	  double x = double(j+1)/double(Nelem)*L;
	  double Qx = Q*(x-L/2);
	  ierr = VecSetValues(b,1,&j,&Qx,INSERT_VALUES); CHKERRA(ierr);
	  // double val = x*(L-x)*Q/cond/2.;
	  double val=Q/cond*x*(L-x)*(x-L/2.)/6.;
	  ierr = VecSetValues(xex,1,&j,&val,INSERT_VALUES); CHKERRA(ierr);
	}
      }
      ierr = VecAssemblyBegin(b); CHKERRA(ierr);
      ierr = VecAssemblyEnd(b); CHKERRA(ierr);

      ierr = VecAssemblyBegin(xex); CHKERRA(ierr);
      ierr = VecAssemblyEnd(xex); CHKERRA(ierr);

      // A.view(VIEWER_STDOUT_WORLD);
      ierr = A.solve(b,x); CHKERRA(ierr); 

      if (debug_print) {
	PetscPrintf(PETSC_COMM_WORLD,"b: \n");
	ierr = VecView(b,VIEWER_STDOUT_WORLD);
	PetscPrintf(PETSC_COMM_WORLD,"FEM: \n");
	ierr = VecView(x,VIEWER_STDOUT_WORLD);
	PetscPrintf(PETSC_COMM_WORLD,"Exact: \n");
	ierr = VecView(xex,VIEWER_STDOUT_WORLD);
      }

      double norm,normex;
      double scal = -1.;
      ierr  = VecNorm(xex,NORM_2,&normex); CHKERRA(ierr);

      ierr = VecAXPY(&scal,xex,x); CHKERRA(ierr);
      ierr  = VecNorm(x,NORM_2,&norm); CHKERRA(ierr);

      PetscPrintf(PETSC_COMM_WORLD,"||x-xex|| = %g,   "
		  "||x-xex||/||xex|| = %g\n",norm,norm/normex);
    }
    A.clean_factor(); 
  }
  PetscFinalize();
  exit(0);
 
}
