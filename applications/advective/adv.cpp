//__INSERT_LICENSE__
//$Id: adv.cpp,v 1.13 2002/09/05 20:25:47 mstorti Exp $
 
#include <src/fem.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/utils.h>
#include <src/util2.h>
#include "advective.h"
#include <time.h>

static char help[] = "Basic finite element program.\n\n";

extern int MY_RANK,SIZE;
int print_internal_loop_conv_g=0,
  consistent_supg_matrix_g=0,
  local_time_step_g=0,
  comp_mat_each_time_step_g=0;

// PETSc now doesn't have the string argument that represents the variable name
// so that I will use this wrapper until I find how to set names in Ascii matlab viewers.
#define PetscViewerSetFormat_WRAPPER(viewer,format,name) \
          PetscViewerSetFormat(viewer,format)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int MyKSPMonitor(KSP ,int ,double ,void *)"
int MyKSPMonitor(KSP ksp,int n,double rnorm,void *dummy) {
  if (print_internal_loop_conv_g) 
    PetscPrintf(PETSC_COMM_WORLD,
		"iteration %d KSP Residual norm %7.4e \n",n,rnorm);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "bless_elemset"
void bless_elemset(char *type,Elemset *& elemset) {
  SET_ELEMSET_TYPE(volume_euler)
    SET_ELEMSET_TYPE(absorb_euler)
    SET_ELEMSET_TYPE(bcconv_adv_euler)
    SET_ELEMSET_TYPE(bcconv_adv_eulerfm2)
    SET_ELEMSET_TYPE(volume_eulerfm2)
    // SET_ELEMSET_TYPE(absorb_eulerfm2)
    SET_ELEMSET_TYPE(volume_shallow)
    SET_ELEMSET_TYPE(absorb_shallow)
    SET_ELEMSET_TYPE(bcconv_adv_shallow)
    SET_ELEMSET_TYPE(bcconv_adv_shallowfm2)
    SET_ELEMSET_TYPE(bcconv_adv_shallowfm2t)
    SET_ELEMSET_TYPE(volume_shallowfm2)
    SET_ELEMSET_TYPE(volume_shallowfm2t)
    {
      printf("not known elemset \"type\": %s\n",type);
      exit(1);
    }
}

#define VECVIEW(name,label) \
ierr = ViewerSetFormat(matlab, \
		       PETSC_VIEWER_ASCII_MATLAB,#label); \
ierr = VecView(name,matlab); CHKERRA(ierr)

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
#undef __FUNC__
#define __FUNC__ "main"
int main(int argc,char **args) {

  Vec     x, dx, xold, res; /* approx solution, RHS, residual*/
  Mat     A_mass;                              /* linear system matrix */
  SLES    sles_mass;     /* linear solver context */
  PC      pc_mass;           /* preconditioner context */
  KSP     ksp_mass;        /* Krylov subspace method context */
  double  norm, *sol, scal; /* norm of solution error */
  int     ierr, i, n = 10, col[3], its, size, node,
    jdof, k, kk, nfixa,
    kdof, ldof, lloc, ndim, nel, nen, neq, nu,
    myrank;
  PetscTruth flg;
  // nu:= dimension of the state vector per node
  PetscScalar  neg_one = -1.0, one = 1.0, value[3];
  PetscScalar *px;
  char fcase[FLEN+1];
  Darray *da; // este me parece que se puede sacar!!
  //Elemset *elemset;
  Dofmap *dofmap;
  dofmap = new Dofmap;
  Mesh *mesh;
  vector<double> dtmin(1,0.);
  Vec a_mass;

  // euler_volume::set_flux_fun(&flux_fun_euler);
  // euler_absorb::flux_fun = &flux_fun_euler;

  // elemsetlist =  da_create(sizeof(Elemset *));
  PetscInitialize(&argc,&args,(char *)0,help);
  print_copyright();

  // Start registering functions
  // Amplitude::initialize_function_table();

  // Get MPI info
  MPI_Comm_size(PETSC_COMM_WORLD,&SIZE);
  MPI_Comm_rank(PETSC_COMM_WORLD,&MY_RANK);

//    MPI_Comm_size(PETSC_COMM_WORLD,&size);
//    MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);

      //  if (size != 1) SETERRA(1,0,"This is a uniprocessor example only!");
  ierr = PetscOptionsGetString(PETSC_NULL,"-case",fcase,FLEN,&flg); CHKERRA(ierr);
  if (!flg) {
    PetscPrintf(PETSC_COMM_WORLD,
		"Option \"-case <filename>\" not passed to PETSc-FEM!!\n");
    PetscFinalize();
    exit(0);
  }

  // Read data
  read_mesh(mesh,fcase,dofmap,neq,SIZE,MY_RANK);
  GLOBAL_OPTIONS = mesh->global_options;

  //o Absolute tolerance when solving a consistent matrix
  GETOPTDEF(double,atol,1e-6);
  //o Relative tolerance when solving a consistent matrix
  GETOPTDEF(double,rtol,1e-3);
  //o Divergence tolerance when solving a consistent matrix
  GETOPTDEF(double,dtol,1e+3);
  //o Maximum iterations when solving a consistent matrix
  GETOPTDEF(int,maxits,150);
  //o Prints the convergence history when solving a consistent matrix
  GETOPTDEF(int,print_internal_loop_conv,0);
  print_internal_loop_conv_g=print_internal_loop_conv;
  //o Measure performance of the \verb+comp_mat_res+ jobinfo. 
  GETOPTDEF(int,measure_performance,0);

  //o Save state vector frequency (in steps)
  GETOPTDEF(int,nsave,10);
  //o Save state vector frequency with the ``rotary save''
  // mechanism. (see \ref{sec:rotary_save})
  GETOPTDEF(int,nsaverot,100);
  //o Sets the number of states saved in a given file
  // in the ``rotary save'' mechanism (see \ref{sec:rotary_save}
  GETOPTDEF(int,nrec,1000000);
  //o Sets the number of files in the ``rotary save'' mechanism. 
  // (see \ref{sec:rotary_save})
  GETOPTDEF(int,nfile,1);

  //o The number of time steps. 
  GETOPTDEF(int,nstep,10000);
  //o Output CPU time statistics for frequency in time steps. 
  GETOPTDEF(int,nstep_cpu_stat,10);
  //o After computing the linear system prints Jacobian and
  // right hand side and stops.. 
  GETOPTDEF(int,print_linear_system_and_stop,0);

  // warning: passed to advective.cpp via a global variable
  //o Uses consistent SUPG matrix for the temporal term or not. 
  GETOPTDEF(int,consistent_supg_matrix,0);
  consistent_supg_matrix_g = consistent_supg_matrix;
  //o Chooses automatically the time step from the 
  // selected Courant number
  GETOPTDEF(int,auto_time_step,1);
  // warning: passed to advective.cpp via a global variable
  //o Chooses a time step that varies locally. (Only makes sense
  // when looking for steady state solutions. 
  GETOPTDEF(int,local_time_step,1);
  local_time_step_g=local_time_step;
  //o The Courant number.
  GETOPTDEF(double,Courant,0.6);
  //o Time step. 
  GETOPTDEF(double,Dt,0.);

  comp_mat_each_time_step_g = 
    consistent_supg_matrix || local_time_step;

  //o Counts time from here.
  GETOPTDEF(double,start_time,0.);
  //o Tolerance when solving with the mass matrix. 
  GETOPTDEF(double,tol_mass,1e-3);

  Time time;
  time.set(start_time);

  //o The pattern to generate the file name to save in for
  // the rotary save mechanism.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file_pattern,outvector%d.out);
  // char *save_file_pattern;
  // mesh->global_options->get_entry("save_file_pattern",save_file_pattern);
  // if (save_file_pattern == NULL) save_file_pattern = "outvector%d.sal";

  //o Filename for saving the state vector.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file,outvector.out);
  // char *save_file;
  // mesh->global_options->get_entry("save_file",save_file);
  // if (save_file == NULL) save_file = "outvector.sal";

#if 0
  PetscViewer matlab;
  ierr = ViewerASCIIOpen(PETSC_COMM_WORLD,
			 "matns.m",&matlab); CHKERRA(ierr);
#endif
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

  dofmap->create_MPI_vector(x);

  // initialize vectors
  ierr = VecDuplicate(x,&xold); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx); CHKERRA(ierr);
  ierr = VecDuplicate(x,&res); CHKERRA(ierr);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // initialize state vectors
  scal=0;
  ierr = VecSet(&scal,x); CHKERRA(ierr);

  arg_list argl;

  if (comp_mat_each_time_step_g) {

    VOID_IT(argl);
    argl.arg_add(&A_mass,PROFILE);
    ierr = assemble(mesh,argl,dofmap,"comp_mat_mass",&time); CHKERRA(ierr);

  } else {

    //DIAG_MAT_MATRIX is obsolete. It was supposed that with this
    //flag the matrix were assembled as a vector and then we used
    // vector_divide(). But this interacts bad with slip boundary
    // conditions. It's better to let the matrix be PETSc matrix
    // (not a vector) and then build the mass matrix ass lumped at
    // the element level. 
#ifdef DIAG_MAT_MATRIX
    ierr = VecDuplicate(x,&a_mass); CHKERRA(ierr);

    VOID_IT(argl);
    argl.arg_add(&a_mass,OUT_VECTOR);
    ierr = assemble(mesh,argl,dofmap,
		    "comp_diag_mat_mass",&time); CHKERRA(ierr);

#else 

    VOID_IT(argl);
    argl.arg_add(&A_mass,PROFILE);
    ierr = assemble(mesh,argl,dofmap,"comp_mat_mass",&time); CHKERRA(ierr);

	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    ierr = SLESCreate(PETSC_COMM_WORLD,&sles_mass); CHKERRA(ierr);
    ierr = SLESSetOperators(sles_mass,A_mass,A_mass,
			    DIFFERENT_NONZERO_PATTERN); CHKERRA(ierr);
    ierr = SLESGetKSP(sles_mass,&ksp_mass); CHKERRA(ierr);
    ierr = SLESGetPC(sles_mass,&pc_mass); CHKERRA(ierr);

    ierr = KSPSetType(ksp_mass,KSPCG); CHKERRA(ierr);
    ierr = PCSetType(pc_mass,PCJACOBI); CHKERRA(ierr);
    ierr = KSPSetTolerances(ksp_mass,tol_mass,PETSC_DEFAULT,PETSC_DEFAULT,
			    PETSC_DEFAULT); CHKERRA(ierr);
    ierr = KSPSetMonitor(ksp_mass,MyKSPMonitor,PETSC_NULL,NULL);

	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    VOID_IT(argl);
    argl.arg_add(&A_mass,OUT_MATRIX);
    ierr = assemble(mesh,argl,dofmap,"comp_mat_mass",&time); CHKERRA(ierr);
#endif
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  ierr = opt_read_vector(mesh,x,dofmap,MY_RANK); CHKERRA(ierr);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

  Chrono chrono; 
#define STAT_STEPS 5
  double cpu[STAT_STEPS],cpuav;
  for (int tstep=1; tstep<=nstep; tstep++) {
    // Take cputime statistics
    if (MY_RANK==0) {
      if (tstep>1) {
	int last=(tstep-2) % STAT_STEPS;
	cpu[last] = chrono.elapsed();
	if (tstep==2) {
	  for (int jstep=1; jstep<STAT_STEPS; jstep++) {
	    cpu[jstep] = cpu[0];
	  }
	}
	cpuav = 0;
	for (int jstep=0; jstep<STAT_STEPS; jstep++) 
	  cpuav+= cpu[jstep];
	cpuav /= STAT_STEPS;
	if (tstep % nstep_cpu_stat == 0)
	  printf("CPU per iter: %f, average [over last %d]: %f\n",
		 cpu[last],STAT_STEPS,cpuav);
      }
    }
    chrono.start();
    ierr = VecCopy(x,xold);

    scal=0;
    ierr = VecSet(&scal,res); CHKERRA(ierr);


    if (comp_mat_each_time_step_g) {

      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
      ierr = SLESCreate(PETSC_COMM_WORLD,&sles_mass); CHKERRA(ierr);
      ierr = SLESSetOperators(sles_mass,A_mass,A_mass,
			      DIFFERENT_NONZERO_PATTERN); CHKERRA(ierr);
      ierr = SLESGetKSP(sles_mass,&ksp_mass); CHKERRA(ierr);
      ierr = SLESGetPC(sles_mass,&pc_mass); CHKERRA(ierr);

      ierr = KSPSetType(ksp_mass,KSPGMRES); CHKERRA(ierr);
      ierr = PCSetType(pc_mass,PCJACOBI); CHKERRA(ierr);
      // ierr = KSPSetTolerances(ksp_mass,tol_mass,PETSC_DEFAULT,PETSC_DEFAULT,
      // PETSC_DEFAULT); CHKERRA(ierr);
      ierr = KSPSetTolerances(ksp_mass,rtol,atol,dtol,maxits); CHKERRA(ierr);
      ierr = KSPSetMonitor(ksp_mass,MyKSPMonitor,PETSC_NULL,NULL); CHKERRA(ierr);
      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

      ierr = MatZeroEntries(A_mass); CHKERRA(ierr);
      VOID_IT(argl);
      argl.arg_add(&x,IN_VECTOR);
      argl.arg_add(&res,OUT_VECTOR);
      argl.arg_add(&dtmin,VECTOR_MIN);
      argl.arg_add(&A_mass,OUT_MATRIX);

#if 0 // In order to measure performance
      Chrono chrono;
      chrono.start();
      int N=10000,nelem = 861;
      for (int jjj=0; jjj<N; jjj++) {
	printf("[loop iter %d]\n",jjj+1);
	ierr = assemble(mesh,argl,dofmap,"comp_res",&time); CHKERRA(ierr);
      }
      double cpu = chrono.elapsed();
      printf("total %f, ntimes %d, nelems %d, rate %f [sec/1000/elems/iter]\n",
	     cpu, N, nelem, cpu/(N*nelem)*1000);
      PetscFinalize();
      exit(0);
#endif 

      if (measure_performance) {
	ierr = measure_performance_fun(mesh,argl,dofmap,"comp_res",
				       &time); CHKERRA(ierr);
	PetscFinalize();
	exit(0);
      }
      ierr = assemble(mesh,argl,dofmap,"comp_res",&time); CHKERRA(ierr);

      ierr = SLESSolve(sles_mass,res,dx,&its); CHKERRA(ierr); 
      ierr = SLESDestroy(sles_mass);

    } else {

      VOID_IT(argl);
      argl.arg_add(&x,IN_VECTOR);
      argl.arg_add(&res,OUT_VECTOR);
      argl.arg_add(&dtmin,VECTOR_MIN);

      if (measure_performance) {
	ierr = measure_performance_fun(mesh,argl,dofmap,"comp_res",
				       &time); CHKERRA(ierr);
	PetscFinalize();
	exit(0);
      }
      ierr = assemble(mesh,argl,dofmap,"comp_res",&time); CHKERRA(ierr);

#ifdef DIAG_MAT_MATRIX
      ierr = VecCopy(res,dx);
      vector_divide(dx,a_mass);
#else
      ierr = SLESSolve(sles_mass,res,dx,&its); CHKERRA(ierr); 
#endif
    }

    // Prints residual and mass matrix in Matlab format
    if (print_linear_system_and_stop) {
      PetscPrintf(PETSC_COMM_WORLD,
		  "Printing residual and matrix for debugging and stopping..\n");
      PetscViewer matlab;
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,
			     "mat.output",&matlab); CHKERRA(ierr);
      ierr = PetscViewerSetFormat_WRAPPER(matlab, 
			     PETSC_VIEWER_ASCII_MATLAB,"res");
      ierr = VecView(res,matlab);
      ierr = PetscViewerSetFormat_WRAPPER(matlab, 
			     PETSC_VIEWER_ASCII_MATLAB,"amass");
      ierr = MatView(A_mass,matlab);
      PetscFinalize();
      exit(0);
    }

    // Define time step depending on strategy. Automatic time step,
    // local time step, etc... 
    if (auto_time_step) Dt = Courant*dtmin[0];
    if (local_time_step) Dt = Courant;
    if (Dt<=0.) {
      PFEMERRQ("error: Dt<=0. You have to set Dt>0 or use "
	       "the auto_time_step flag.\n");
    }
    // SHV(Dt);
    time.inc(Dt);

    // Upgrade state vector.
    ierr = VecAXPY(&Dt,dx,x);

    VOID_IT(argl);
    argl.arg_add(&x,IN_OUT_VECTOR);
    argl.arg_add(&xold,IN_VECTOR);
    ierr = assemble(mesh,argl,dofmap,"absorb_bc_proj",&time); CHKERRA(ierr);

    ierr  = VecNorm(res,NORM_2,&norm); CHKERRA(ierr);
    double time_=time;

    // Reuse dx for computing the delta state vector after re-projecting
    double delta_u;
    ierr = VecCopy(x,dx);
    scal=-1.;
    ierr = VecAXPY(&scal,xold,dx);
    ierr  = VecNorm(dx,NORM_2,&delta_u); CHKERRA(ierr);

    PetscPrintf(PETSC_COMM_WORLD,
		"time_step %d, time: %g, res = %14.12e - delta_u = %10.3e\n",
		tstep,time_,norm,delta_u);
    print_vector_rota(save_file_pattern.c_str(),x,dofmap,
		      &time,tstep-1,nsaverot,nrec,nfile);

    if (tstep % nsave == 0){
      PetscPrintf(PETSC_COMM_WORLD,
		  " --------------------------------------\n"
		  "Time step: %d\n"
		  " --------------------------------------\n",
		  tstep);

      print_vector(save_file.c_str(),x,dofmap,&time);
    }
  }
  print_vector(save_file.c_str(),x,dofmap,&time);

  ierr = VecDestroy(x); CHKERRA(ierr); 
  ierr = VecDestroy(xold); CHKERRA(ierr); 
  ierr = VecDestroy(dx); CHKERRA(ierr); 
  ierr = VecDestroy(res); CHKERRA(ierr); 
#ifdef DIAG_MAT_MATRIX
  ierr = MatDestroy(A_mass); CHKERRA(ierr); 
#endif
  
  PetscFinalize();
  exit(0);
}
