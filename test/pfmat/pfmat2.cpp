/*__INSERT_LICENSE__*/
// $Id: pfmat2.cpp,v 1.4 2002/07/29 20:17:39 mstorti Exp $

// Tests for the `PFMat' class

#include <src/debug.h>

#include <stdlib.h>
#include <src/fem.h>
#include <src/utils.h>
#include <src/elemset.h>
#include <src/pfmat.h>
#include <src/pfptscmat.h>
#include <src/iisdmat.h>
//#include <src/petscmat.h>
#include <src/graph.h>

// Runs a simple example for testing the PFMat matrix classes.
// 

double drand_r(unsigned int *SEED) {
  return double(rand_r(SEED))/double(RAND_MAX);
}

extern int MY_RANK,SIZE;

class Part : public DofPartitioner {
public:
  int N,comm_size;
  int processor(int k) const { return k*comm_size/N; }
  ~Part() {}
} part;

#define CHKOPT(name)  if (myrank==0) assert(!strcmp(#name,args[++arg]))
//  void set_print_maybe(int cond, int myrank, PFMat &A,int j, int k,
//  		     double val) {
  
#define SET_PRINT_MAYBE(j,k,val)					\
     { if (((j)*size)/Nelem == myrank) A.set_value((j),(k),(val));	\
       if (myrank == 0) fprintf(fid,"%d %d %f\n",(j),(k),(val)); }

int main(int argc,char **args) {

  FILE *fid;
  PFMat *A_p;
  Vec x,b,xex;
  int ierr,iisd_subpart,nsolve,rand_flag,mat_type,q_type,
    tests_ok=1,nprof=1,chkkey,chkkey_c,oct_check=1;
  unsigned int SEED;
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
  // <iisd_subpart> nmat <nmat> mtyp <mat_type> f <fill>

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

  CHKOPT(nprof);
  nprof=1;
  if (myrank==0 && argc>++arg) 
    sscanf(args[arg],"%d",&nprof);
  ierr = MPI_Bcast (&nprof, 
		    1, MPI_INT, 0,MPI_COMM_WORLD); CHKERRA(ierr); 
  PetscPrintf(PETSC_COMM_WORLD,"nprof %d  (solve with nprof profiles)\n",nmat);

  rand_flag=0;

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
    sscanf(args[arg],"%lf",&fill);
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
  //  PETScMat AAA(N,N,part,PETSC_COMM_WORLD);
  if (mat_type==0) {
    A_p = &AA;
  } else if (mat_type==1) {
    // A_p = &AAA;
    assert(0);
  } else assert(0);
  PFMat &A = *A_p;

  A.set_option("iisd_subpart",iisd_subpart);
  A.set_option("print_internal_loop_conv",debug_print);
  A.set_option("rtol",1e-8);
  A.set_option("atol",0);
  A.set_option("print_fsm_transition_info",0);
  A.set_option("use_compact_profile",0);

  for (int jprof=0; jprof<nprof; jprof++) {

    // define type of elements
    is_diag.clear();
    is_diag.resize(Nelem,0);
    if (myrank==0) {
      for (int je=0; je<Nelem; je++) {
	double v = drand_r(&SEED);
	// printf("%d %f\n",je,v);
	is_diag[je] = v > fill;
      }
    }
    ierr = MPI_Bcast (is_diag.begin(), Nelem, MPI_INT, 0,MPI_COMM_WORLD);
    
    // Define profile
    for (int j=0; j<Nelem; j++) {
      // Is this element in this processor?
      if ((j*size)/Nelem != myrank) continue; 
      if (j-1>=0) A.set_profile(j-1,j-1);
      if (j<N)   A.set_profile(j,j);
      if (!is_diag[j] && j-1>=0 && j<N) {
	A.set_profile(j-1,j);
	A.set_profile(j,j-1);
      }
    }

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
      // A.zero_entries();

      if (myrank==0) {
	fid = fopen("data.tmp","w");
	chkkey = rand();
	fprintf(fid,"%d %d %g\n",chkkey,0,0.);
      }
	
      for (int j=0; j<Nelem; j++) {
	// Is this element in this processor?
	if (j-1>=0) SET_PRINT_MAYBE(j-1,j-1,1.);
	if (j<N)   SET_PRINT_MAYBE(j,j,1.);
	if (!is_diag[j] && j-1>=0 && j<N) {
	  SET_PRINT_MAYBE(j-1,j,-1.);
	  SET_PRINT_MAYBE(j,j-1,-1.);
	}
      }
      A.assembly_begin(MAT_FINAL_ASSEMBLY); 
      A.assembly_end(MAT_FINAL_ASSEMBLY); 
      if (myrank==0) {
	fclose(fid);
	if (oct_check) system("octave -q < checkpf.m >/dev/null");
      }
      // A.view();

      for (int ksolve=0; ksolve<nsolve; ksolve++) {

	if (!oct_check) PetscPrintf(PETSC_COMM_WORLD,"solve...\n");
	if (myrank==0) {
	  for (int j=0; j<N; j++) {
	    double Qx=1.;
	    ierr = VecSetValues(b,1,&j,&Qx,INSERT_VALUES);
	    CHKERRA(ierr);
	  }
	}
	if (oct_check && myrank==0) {
	  fid = fopen("xsol.tmp","r");
	  fscanf(fid,"%d",&chkkey_c);
	  assert(chkkey_c==chkkey);
	  for (int j=0; j<N; j++) {
	    double val;
	    fscanf(fid,"%lf",&val);
	    ierr = VecSetValues(xex,1,&j,&val,INSERT_VALUES);
	    CHKERRA(ierr);
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
	  if (oct_check) {
	    PetscPrintf(PETSC_COMM_WORLD,"Exact (Octave): \n");
	    ierr = VecView(xex,VIEWER_STDOUT_WORLD);
	  }
	}

	if (oct_check) {
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
      }
      A.clean_factor(); 
    }
    A.clear();
  }

  if (oct_check) PetscPrintf(PETSC_COMM_WORLD,
			     "All tests OK ?  %d\n",tests_ok);
  PetscFinalize();
}
