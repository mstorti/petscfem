/*
  This file belongs to he PETSc - FEM package a library and
  application suite oriented to the Finite Element Method based on PETSc. 
  Copyright (C) 1999, 2000  Mario Alberto Storti
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.

*/

#include <set>

#include "../../src/fem.h"
#include "../../src/readmesh.h"
#include "../../src/getprop.h"
#include "../../src/utils.h"
#include "../../src/util2.h"
#include "advective.h"
#include <time.h>

static char help[] = "Basic finite element program.\n\n";

int MY_RANK,SIZE;
TextHashTable *GLOBAL_OPTIONS;
int print_internal_loop_conv_g=0,
  consistent_supg_matrix_g=0,
  local_time_step_g=0,
  comp_mat_each_time_step_g=0;

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
#if 0
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
    SET_ELEMSET_TYPE(volume_shallowfm2)
    SET_ELEMSET_TYPE(volume_shallowfm2t)
#endif
    //    SET_ELEMSET_TYPE(bcconv_adv_advecfm2)
    // SET_ELEMSET_TYPE(advdif_advecfm2)
    if ( ! strcmp(type,"advdif_advecfm2")) { 
      advdif_advecfm2 * pp = new advdif_advecfm2; 
      elemset = pp;
    } else
//      SET_ELEMSET_TYPE(bcconv_adv_burgers)
//      SET_ELEMSET_TYPE(advdif_burgers)
    {
      printf("not known elemset \"type\": %s\n",type);
      exit(1);
    }
}

#define VECVIEW(name,label) \
ierr = ViewerSetFormat(matlab, \
		       VIEWER_FORMAT_ASCII_MATLAB,#label); \
ierr = VecView(name,matlab); CHKERRA(ierr)

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
#undef __FUNC__
#define __FUNC__ "main"
int main(int argc,char **args) {

  Vec     x, dx, xold, res; /* approx solution, RHS, residual*/
  Mat     A;		/* linear system matrix */
#ifdef CHECK_JAC
  Mat     AA;
#endif
  SLES    sles;		/* linear solver context */
  PC      pc;		/* preconditioner context */
  KSP     ksp;		/* Krylov subspace method context */
  double  *sol, scal;	/* norm of solution error */
  int     ierr, i, n = 10, col[3], its, flg, size, node,
    jdof, k, kk, nfixa,
    kdof, ldof, lloc, ndim, nel, nen, neq, nu,
    myrank;
  // nu:= dimension of the state vector per node
  Scalar  neg_one = -1.0, one = 1.0, value[3];
  Scalar *px;
  char fcase[FLEN+1];
  Darray *da; // este me parece que se puede sacar!!
  //Elemset *elemset;
  Dofmap *dofmap;
  dofmap = new Dofmap;
  Mesh *mesh;
  vector<double> dtmin(1,0.);
  Vec a;
  GlobParam glob_param;

  // euler_volume::set_flux_fun(&flux_fun_euler);
  // euler_absorb::flux_fun = &flux_fun_euler;

  // elemsetlist =  da_create(sizeof(Elemset *));
  PetscInitialize(&argc,&args,(char *)0,help);
  print_copyright();

  // Start registering functions
  Amplitude::initialize_function_table();

  // Get MPI info
  MPI_Comm_size(PETSC_COMM_WORLD,&SIZE);
  MPI_Comm_rank(PETSC_COMM_WORLD,&MY_RANK);

//    MPI_Comm_size(PETSC_COMM_WORLD,&size);
//    MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);

      //  if (size != 1) SETERRA(1,0,"This is a uniprocessor example only!");
  ierr = OptionsGetString(PETSC_NULL,"-case",fcase,FLEN,&flg); CHKERRA(ierr);

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
  //o Solve system before \verb+print\_linear_system_and_stop+
  GETOPTDEF(int,solve_system,1);

  //o Sets the save frequency in iterations for the ``print some''
  // mechanism. (see doc in the Navier-Stokes module)
  GETOPTDEF(int,nsome,10000);
  //o Name of file where to read the nodes for the ``print some'' 
  // feature. 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,print_some_file,);
  //o Name of file where to save node values for the ``print some'' 
  // feature. 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file_some,outvsome.out);

  set<int> node_list;
  print_some_file_init(mesh->global_options,
		       print_some_file.c_str(),
		       save_file_some.c_str(),node_list);


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
  glob_param.Dt = Dt;
  //o The parameter of the trapezoidal rule
  // for temporal integration. 
  GETOPTDEF(double,alpha,0.);
  glob_param.alpha=alpha;
#define ALPHA (glob_param.alpha)
  if (ALPHA>0.) {
    consistent_supg_matrix = 1;
    auto_time_step=0;		// Don't know how to do auto_time_step
				// in the implicit case...
    local_time_step=0;
  }
  //o Number of iterations in the Newton loop. (
  // for the implicit method. 
  GETOPTDEF(int,nnwt,3);
  if (ALPHA==0.) nnwt=1;

  comp_mat_each_time_step_g = 
    consistent_supg_matrix || local_time_step;

  //o Counts time from here.
  GETOPTDEF(double,start_time,0.);
  //o Tolerance when solving with the mass matrix. 
  GETOPTDEF(double,tol_mass,1e-3);
  //o Tolerance when solving the sublinear problem
  // at each iteration.
  GETOPTDEF(double,tol_linear,1e-3);
  //o Tolerance when solving the non-linear problem
  // for the implicit case.
  GETOPTDEF(double,tol_newton,1e-3);

  //o Chooses the preconditioning operator. 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,preco_type,Jacobi);
  // I had to do this since `c_str()' returns `const char *'
  char *preco_type_ = new char[preco_type.size()+1];
  strcpy(preco_type_,preco_type.c_str());

  Time time,time_star;
  time.set(start_time);

  //o The pattern to generate the file name to save in for
  // the rotary save mechanism.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file_pattern,outvector%d.out);

  //o Filename for saving the state vector.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file,outvector.out);

#if 0
  Viewer matlab;
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

  VOID_IT(argl);
  argl.arg_add(&A,PROFILE);
  ierr = assemble(mesh,argl,dofmap,"comp_prof",&time); CHKERRA(ierr);

#ifdef CHECK_JAC
  VOID_IT(argl);
  argl.arg_add(&AA,PROFILE);
  ierr = assemble(mesh,argl,dofmap,"comp_prof",&time); CHKERRA(ierr);
#endif

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  ierr = opt_read_vector(mesh,x,dofmap,MY_RANK); CHKERRA(ierr);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // This is for taking statistics of the
  // CPU time consumed by a time steptime
  Chrono chrono; 
#define STAT_STEPS 5
  double cpu[STAT_STEPS],cpuav;
  for (int tstep=1; tstep<=nstep; tstep++) {
    time_star.set(time.time()+alpha*Dt);
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

    for (int inwt=0; inwt<nnwt; inwt++) {

      // Initializes res
      scal=0;
      ierr = VecSet(&scal,res); CHKERRA(ierr);

      if (comp_mat_each_time_step_g) {

	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
	ierr = SLESCreate(PETSC_COMM_WORLD,&sles); CHKERRA(ierr);
	ierr = SLESSetOperators(sles,A,A,
				DIFFERENT_NONZERO_PATTERN); CHKERRA(ierr);
	ierr = SLESGetKSP(sles,&ksp); CHKERRA(ierr);
	ierr = SLESGetPC(sles,&pc); CHKERRA(ierr);
	
	ierr = KSPSetType(ksp,KSPGMRES); CHKERRA(ierr);
//  	ierr = PCSetType(pc,
//  			 (preco_type=="Jacobi" ? PCJACOBI :
//  			  preco_type=="none" ? PCNONE : 
//  			  preco_type=="LU" ? PCLU : 
//  			  PCJACOBI)); CHKERRA(ierr);
	ierr = KSPSetPreconditionerSide(ksp,PC_RIGHT);
	ierr = PCSetType(pc,preco_type_); CHKERRA(ierr);
	// ierr = PCSetType(pc,PCJACOBI); CHKERRA(ierr);
	ierr = KSPSetTolerances(ksp,rtol,atol,dtol,maxits); CHKERRA(ierr);
	ierr = KSPSetMonitor(ksp,MyKSPMonitor,PETSC_NULL); CHKERRA(ierr);
	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

	ierr = MatZeroEntries(A); CHKERRA(ierr);
#ifdef CHECK_JAC
	ierr = MatZeroEntries(AA); CHKERRA(ierr);
#endif
	VOID_IT(argl);
	argl.arg_add(&xold,IN_VECTOR);
#ifndef CHECK_JAC
	argl.arg_add(&x,IN_VECTOR);
#else
	argl.arg_add(&x,PERT_VECTOR);
#endif
	argl.arg_add(&res,OUT_VECTOR);
	argl.arg_add(&dtmin,VECTOR_MIN);
	argl.arg_add(&A,OUT_MATRIX);
	argl.arg_add(&glob_param,USER_DATA);
#ifdef CHECK_JAC
	argl.arg_add(&AA,OUT_MATRIX_FDJ);
#endif

	if (measure_performance) {
	  ierr = measure_performance_fun(mesh,argl,dofmap,"comp_res",
					 &time_star); CHKERRA(ierr);
	  PetscFinalize();
	  exit(0);
	}
	ierr = assemble(mesh,argl,dofmap,"comp_res",&time_star); CHKERRA(ierr);

	if (!print_linear_system_and_stop || solve_system) {
	  ierr = SLESSolve(sles,res,dx,&its); CHKERRA(ierr); 
	}
	ierr = SLESDestroy(sles);

      } else {

	VOID_IT(argl);
	argl.arg_add(&x,IN_VECTOR);
	argl.arg_add(&res,OUT_VECTOR);
	argl.arg_add(&dtmin,VECTOR_MIN);

	if (measure_performance) {
	  ierr = measure_performance_fun(mesh,argl,dofmap,"comp_res",
					 &time_star); CHKERRA(ierr);
	  PetscFinalize();
	  exit(0);
	}
	ierr = assemble(mesh,argl,dofmap,"comp_res",&time_star); CHKERRA(ierr);

	if (!print_linear_system_and_stop || solve_system) {
	  ierr = SLESSolve(sles,res,dx,&its); CHKERRA(ierr); 
	}
      }

      if (print_linear_system_and_stop) {
	PetscPrintf(PETSC_COMM_WORLD,
		    "Printing residual and matrix for debugging and stopping..\n");
	Viewer matlab;
	ierr = ViewerASCIIOpen(PETSC_COMM_WORLD,
			       "mat.output",&matlab); CHKERRA(ierr);
	ierr = ViewerSetFormat(matlab, 
			       VIEWER_FORMAT_ASCII_MATLAB,"res");
	ierr = VecView(res,matlab);
	if (solve_system) {
	  ierr = ViewerSetFormat(matlab, 
				 VIEWER_FORMAT_ASCII_MATLAB,"dx");
	  ierr = VecView(dx,matlab);
	}
	ierr = ViewerSetFormat(matlab, 
			       VIEWER_FORMAT_ASCII_MATLAB,"A");
	ierr = MatView(A,matlab);
#ifdef CHECK_JAC
	ierr = ViewerSetFormat(matlab, 
			       VIEWER_FORMAT_ASCII_MATLAB,"AA");
	ierr = MatView(AA,matlab);
#endif
	PetscFinalize();
	exit(0);
      }

      double normres;
      ierr  = VecNorm(res,NORM_2,&normres); CHKERRA(ierr);
      PetscPrintf(PETSC_COMM_WORLD,
		  "Newton subiter %d, norm_res  = %10.3e\n",
		  inwt,normres);
      scal=1.;
      ierr = VecAXPY(&scal,dx,x);
      if (normres < tol_newton) break;
    }

    // Prints residual and mass matrix in Matlab format
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

    VOID_IT(argl);
    argl.arg_add(&x,IN_OUT_VECTOR);
    argl.arg_add(&xold,IN_VECTOR);
    ierr = assemble(mesh,argl,dofmap,"absorb_bc_proj",&time_star); CHKERRA(ierr);

    double time_=time;

    // Reuse dx for computing the delta state vector after re-projecting
    double delta_u;
    ierr = VecCopy(x,dx);
    scal=-1.;
    ierr = VecAXPY(&scal,xold,dx);
    ierr  = VecNorm(dx,NORM_2,&delta_u); CHKERRA(ierr);

    PetscPrintf(PETSC_COMM_WORLD,
		"time_step %d, time: %g, delta_u = %10.3e\n",
		tstep,time_,delta_u);
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
    if (print_some_file!="" && tstep % nsome == 0)
      print_some(save_file_some.c_str(),x,dofmap,node_list,&time);

  }
  print_vector(save_file.c_str(),x,dofmap,&time);

  ierr = VecDestroy(x); CHKERRA(ierr); 
  ierr = VecDestroy(xold); CHKERRA(ierr); 
  ierr = VecDestroy(dx); CHKERRA(ierr); 
  ierr = VecDestroy(res); CHKERRA(ierr); 
#ifdef DIAG_MAT_MATRIX
  ierr = MatDestroy(A); CHKERRA(ierr); 
#endif
  
  PetscFinalize();
  exit(0);
}
