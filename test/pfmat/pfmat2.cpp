/*__INSERT_LICENSE__*/
// $Id: pfmat2.cpp,v 1.1.2.1 2002/01/02 23:55:24 mstorti Exp $

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

// Runs a simple example for testing the PFMat matrix classes.
// 

extern int MY_RANK,SIZE;

class Part : public DofPartitioner {
public:
  int N,comm_size;
  int processor(int k) const { return k*comm_size/N; }
  ~Part() {}
} part;

#define CHKOPT(name)  if (myrank==0) assert(!strcmp(#name,args[++arg]))

int main(int argc,char **args) {

  PFMat *A_p;
  Vec x,b,xex;
  int ierr,iisd_subpart,nsolve,rand_flag,mat_type,q_type,
    tests_ok=1;
  // is_diag:= `is_diag[k]' flags whether the element `k' is Laplace
  // (element matrix propto. [1 -1; -1 1]) or Identity ([1 0; 0 1]);
  vector<int> is_diag;
  double tol=1e-10,fill=1.;
  // debug.activate();
  int myrank,size;
  PetscInitialize(&argc,&args,NULL,NULL);

  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);

  double cond,L,Q,rand_coef=1.;
  int Nelem=10, debug_print=0, nmat;
  int arg=0;

  // usage: pfmat.bin ne <Nelem> deb <debug_print> nsl <nsolve> sbdp
  // <iisd_subpart> nmat <nmat> rrhs <rand_flag> mtyp <mat_type>

  CHKOPT(ne);
  if (myrank==0 && argc>++arg)
    sscanf(args[arg],"%d",&Nelem);
  ierr = MPI_Bcast (&Nelem, 
		    1, MPI_INT, 0,MPI_COMM_WORLD); CHKERRA(ierr); 
  int N=Nelem-1;
  PetscPrintf(PETSC_COMM_WORLD,"Nelem %d  (Number of elements)\n",Nelem);

  CHKOPT(deb);
  if (myrank==0 && argc>++arg)
    sscanf(args[arg],"%d",&debug_print);
  ierr = MPI_Bcast (&debug_print, 
		    1, MPI_INT, 0,MPI_COMM_WORLD); CHKERRA(ierr); 
  PetscPrintf(PETSC_COMM_WORLD,"debug_print %d  "
	      "(print debugging information)\n",debug_print);

  CHKOPT(nsl);
  nsolve=1;
  if (myrank==0 && argc>++arg)
    sscanf(args[arg],"%d",&nsolve);
  ierr = MPI_Bcast (&nsolve, 1, 
		    MPI_INT, 0,MPI_COMM_WORLD); CHKERRA(ierr); 
  PetscPrintf(PETSC_COMM_WORLD,
	      "nsolve %d  (repeats solution stage nsolve times\n",
	      nsolve);
    
  CHKOPT(sbdp);
  iisd_subpart=1;
  if (myrank==0 && argc>++arg) 
    sscanf(args[arg],"%d",&iisd_subpart);
  ierr = MPI_Bcast (&iisd_subpart, 
		    1, MPI_INT, 0,MPI_COMM_WORLD); CHKERRA(ierr); 
  PetscPrintf(PETSC_COMM_WORLD,"iisd_subpart %d  "
	      "(number of partitions for subpartitioning"
	      " inside each processor)\n",iisd_subpart);

  CHKOPT(nmat);
  nmat=1;
  if (myrank==0 && argc>++arg) 
    sscanf(args[arg],"%d",&nmat);
  ierr = MPI_Bcast (&nmat, 
		    1, MPI_INT, 0,MPI_COMM_WORLD); CHKERRA(ierr); 
  PetscPrintf(PETSC_COMM_WORLD,"nmat %d  (solve with nmat matrices)\n",nmat);

  CHKOPT(rrhs);
  rand_flag=0;
  if (myrank==0 && argc>++arg) 
    sscanf(args[arg],"%d",&rand_flag);
  ierr = MPI_Bcast (&rand_flag, 
		    1, MPI_INT, 0,MPI_COMM_WORLD); CHKERRA(ierr); 
  PetscPrintf(PETSC_COMM_WORLD,"rand_flag %d  (use a randomly "
	      "perturbated Q, k and L\n",rand_flag);

  CHKOPT(mtyp);
  mat_type=0; // 0 -> IISDMat, 1 -> PETScMat
  if (myrank==0 && argc>++arg) 
    sscanf(args[arg],"%d",&mat_type);
  ierr = MPI_Bcast (&mat_type, 
		    1, MPI_INT, 0,MPI_COMM_WORLD); CHKERRA(ierr); 
  PetscPrintf(PETSC_COMM_WORLD,
	      "mat_type %d (0 -> IISDMat, 1 -> PETScMat)\n",mat_type);

  CHKOPT(f);
  mat_type=0; // 0 -> IISDMat, 1 -> PETScMat
  if (myrank==0 && argc>++arg) 
    sscanf(args[arg],"%f",&fill);
  ierr = MPI_Bcast (&fill, 
		    1, MPI_DOUBLE, 0,MPI_COMM_WORLD); CHKERRA(ierr); 
  PetscPrintf(PETSC_COMM_WORLD,
	      "fill %f (1 -> Laplace, 0 -> Ident.)\n",fill);

  q_type=0; // Anyway, we have no solution unless `fill=1.'

  part.comm_size = size;
  part.N = N;
  MY_RANK = myrank;
  SIZE = size;

  IISDMat AA(N,N,part,PETSC_COMM_WORLD);
  PETScMat AAA(N,N,part,PETSC_COMM_WORLD);
  if (mat_type==0) {
    A_p = &AA;
  } else if (mat_type==1) {
    A_p = &AAA;
  } else assert(0);
  PFMat &A = *A_p;

  A.set_option("iisd_subpart",iisd_subpart);
  A.set_option("print_internal_loop_conv",debug_print);
  A.set_option("rtol",1e-8);
  A.set_option("atol",0);
  A.create();
  // if (debug_print) A.view();

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
	  double Qx = (q_type==0 ? Q : Q*(x-L/2));
	  ierr = VecSetValues(b,1,&j,&Qx,INSERT_VALUES); CHKERRA(ierr);
	  double val = (q_type==0 ? x*(L-x)*Q/cond/2.
			: Q/cond*x*(L-x)*(x-L/2.)/6.);
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

      int this_ok = norm/normex<=tol;
      PetscPrintf(PETSC_COMM_WORLD,"test OK: %d, ||x-xex|| = %g,   "
		  "||x-xex||/||xex|| = %g, tol %g\n",this_ok,
		  norm,norm/normex,tol);
      if (!this_ok) tests_ok = 0;
    }
    A.clean_factor(); 
  }
  PetscPrintf(PETSC_COMM_WORLD,"All tests OK ?  %d\n",tests_ok);
  PetscFinalize();
  exit(0);
 
}
