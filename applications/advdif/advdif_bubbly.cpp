//__INSERT_LICENSE__
//$Id: advdif_bubbly.cpp,v 1.14.10.1 2007/02/19 20:23:56 mstorti Exp $

#include <src/debug.h>
#include <set>

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/pfmat.h>
#include <src/hook.h>

#include "advective.h"

#include <time.h>

static char help[] = "Basic finite element program.\n\n";

extern int print_internal_loop_conv_g,
  consistent_supg_matrix_g,
  local_time_step_g,
  comp_mat_each_time_step_g;

extern const char * jobinfo_fields;

#define VECVIEW(name,label) \
ierr = PetscViewerSetFormat(matlab, \
		       PETSC_VIEWER_ASCII_MATLAB,#label); \
ierr = VecView(name,matlab); CHKERRA(ierr)

// PETSc now doesn't have the string argument that represents the variable name
// so that I will use this wrapper until I find how to set names in Ascii matlab viewers.
#define PetscViewerSetFormat_WRAPPER(viewer,format,name) \
          PetscViewerSetFormat(viewer,format)

Hook *advdif_hook_factory(const char *name);

//-------<*>-------<*>-------<*>-------<*>-------<*>-------
#undef __FUNC__
#define __FUNC__ "main"
int bubbly_main() {

  Vec     x, dx, xold, res; /* approx solution, RHS, residual*/
  Vec     dx_out,res_out;
  // Vec     dx_liq,res_liq;
//  PFMat *A_liq,*A_gas,*AA_liq,*AA_gas;	   // linear system matrix
//  PFMat *A_stage[2,NULL],*AA_stage[2,NULL];	   // linear system matrix
//  PFMat *A_stage,*AA_stage;	   // linear system matrix
  vector<PFMat *> A_stage,AA_stage;

  double  *sol, scal, normres, normres_ext=NAN;    /* norm of solution error */
  int     ierr, i, n = 10, col[3], its, size, node,
    jdof, k, kk, nfixa,
    kdof, ldof, lloc, ndim, nel, nen, neq, nu,
    myrank;
//  int nnwt_in,nstage=2;
  int nnwt_in=0;
  double omega_newton_in=NAN;


  PetscBool flg;
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
  Vec a;
  GlobParam glob_param;
  GLOB_PARAM = &glob_param;
  string save_file_res;

  // euler_volume::set_flux_fun(&flux_fun_euler);
  // euler_absorb::flux_fun = &flux_fun_euler;

  // elemsetlist =  da_create(sizeof(Elemset *));

  print_copyright();
  PetscPrintf(PETSCFEM_COMM_WORLD,
	      "-------- Generic Advective-Diffusive / Bubbly module ---------\n");

  Debug debug2(0,PETSCFEM_COMM_WORLD);

  ierr = PetscOptionsGetString(PETSC_NULL,"-case",fcase,FLEN,&flg); CHKERRA(ierr);
  if (!flg) {
    PetscPrintf(PETSCFEM_COMM_WORLD,
		"Option \"-case <filename>\" not passed to PETSc-FEM!!\n");
    PetscFinalize();
    exit(0);
  }

  // Read data
  ierr = read_mesh(mesh,fcase,dofmap,neq,SIZE,MY_RANK); CHKERRA(ierr);
  GLOBAL_OPTIONS = mesh->global_options;

  //o number of stages in which the continuum problem is split
  GETOPTDEF(int,nstage,1);

  A_stage.resize(nstage,NULL);
  AA_stage.resize(nstage,NULL);

  //o Activate debugging
  GETOPTDEF(int,activate_debug,0);
  if (activate_debug) {
    debug2.activate();
    Debug::init();
  }
  //o Activate printing in debugging
  GETOPTDEF(int,activate_debug_print,0);
  if (activate_debug_print) debug2.activate("print");
  //o Activate report of memory usage
  GETOPTDEF(int,activate_debug_memory_usage,0);
  if (activate_debug_memory_usage) debug2.activate("memory_usage");

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
  //o Measure performance of the #comp_mat_res# jobinfo.
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
  //o Solve system before #print_linear_system_and_stop#
  GETOPTDEF(int,solve_system,1);
  //o If #print_linear_system_and_stop# is active,
  // then print system in this Newton iteration
  GETOPTDEF(int,inwt_stop,0);
  //o If #print_linear_system_and_stop# is active,
  // then print system in this time step
  GETOPTDEF(int,time_step_stop,1);
  //o Print the residual each #nsave# steps.
  GETOPTDEF(int,print_residual,0);

  //o Sets the save frequency in iterations for the ``print some''
  // mechanism. (see doc in the Navier-Stokes module)
  GETOPTDEF(int,nsome,10000);
  //o Name of file where to read the nodes for the ``print some''
  // feature.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,print_some_file,);
  //o Name of file where to save node values for the ``print some''
  // feature.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file_some,outvsome.out);
  //o Access mode to the ``some'' file. If 0 rewind file. If 1
  //  append to previous  results.
  TGETOPTDEF(GLOBAL_OPTIONS,int,save_file_some_append,1);
  //o Print, after execution, a report of the times a given option
  // was accessed. Useful for detecting if an option was used or not.
  GETOPTDEF(int,report_option_access,1);
  //o Use IISD (Interface Iterative Subdomain Direct) or not.
  GETOPTDEF(int,use_iisd,0);
  //o Type of solver. May be #iisd# or #petsc#
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,solver,petsc);
  if (use_iisd) solver = string("iisd");

  // Use IISD (Interface Iterative Subdomain Direct) or not.
  // A_tet = (use_iisd ? &IISD_A_tet : &PETSc_A_tet);

//  A_liq  = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());
//  AA_liq = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());
//  A_gas  = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());
//  AA_gas = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());

    for (int kstage=0; kstage<nstage; kstage++) {
  A_stage[kstage]  = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());
  AA_stage[kstage] = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());
}

  set<int> node_list;
  print_some_file_init(mesh->global_options,
		       print_some_file.c_str(),
		       save_file_some.c_str(),node_list,
		       save_file_some_append);

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
  // The value of #Dt# is used for printing etc... If #Dt#
  // is not set and #steady+ is set then #Dt# is set to one.
  GETOPTDEF(int,steady,0);
  if (steady && Dt==0.) Dt=1.;
  glob_param.Dt = Dt;
  glob_param.steady = steady;
  //o The parameter of the trapezoidal rule
  // for temporal integration.
  GETOPTDEF(double,alpha,0.);
  glob_param.alpha=alpha;
#define ALPHA (glob_param.alpha)

  //o Flag if steady solution or not (uses Dt=inf). If #steady#
  GETOPTDEF(int,nnwt_liq,3);
  GETOPTDEF(int,nnwt_gas,3);
  // is set to 1, then the computations are as if $\Dt=\infty$.
  if (ALPHA==0.) nnwt_liq=1;
  if (ALPHA==0.) nnwt_gas=1;

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
  //o Tolerance when solving for a steady state
  GETOPTDEF(double,tol_steady,0.);
  //o Relaxation factor for the Newton iteration
  GETOPTDEF(double,omega_newton,1.);
  GETOPTDEF(double,omega_newton_liq,omega_newton);
  GETOPTDEF(double,omega_newton_gas,omega_newton);

  vector<double> gather_values;
  //o Number of ``gathered'' quantities.
  GETOPTDEF(int,ngather,0);
  //o Print values in this file 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,gather_file,gather.out);
  // Initialize gather_file
  FILE *gather_file_f;
  if (MY_RANK==0 && ngather>0) {
    gather_file_f = fopen(gather_file.c_str(),"w");
    //fprintf(gather_file_f,"");
    fclose(gather_file_f);
  }

  //o Chooses the preconditioning operator.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,preco_type,jacobi);
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
  save_file_res = save_file + string(".res");

#if 0
  PetscViewer matlab;
  ierr = PetscViewerASCIIOpen(PETSCFEM_COMM_WORLD,
			      "matns.m",&matlab); CHKERRA(ierr);
#endif
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:

  dofmap->create_MPI_vector(x);

  // initialize vectors
  ierr = VecDuplicate(x,&xold); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx); CHKERRA(ierr);
  ierr = VecDuplicate(x,&res); CHKERRA(ierr);
  ierr = VecDuplicate(x,&res_out); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx_out); CHKERRA(ierr);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // initialize state vectors
  scal=0;
  ierr = VecSet(x,scal); CHKERRA(ierr);

  arg_list argl,arglf;

  for (int kstage=0; kstage<nstage; kstage++) {

    if (nstage==1) {
      jobinfo_fields = "gasliq";
      nnwt_in = nnwt_liq;
      omega_newton_in = omega_newton_liq;
    }

    else if (kstage==0) {
      jobinfo_fields = "liq";
      nnwt_in = nnwt_liq;
      omega_newton_in = omega_newton_liq;
    }
    else if (kstage==1) {
      jobinfo_fields = "gas";
      nnwt_in = nnwt_gas;
      omega_newton_in = omega_newton_gas;
    }
    else if (kstage==2) {
      jobinfo_fields = "kep";
      nnwt_in = nnwt_liq;
      omega_newton_in = omega_newton_liq;
    }

    VOID_IT(argl);
    //  argl.arg_add(A_liq,PROFILE|PFMAT);
    //  jobinfo_fields = "liq";
    argl.arg_add(A_stage[kstage],PROFILE|PFMAT);
    ierr = assemble(mesh,argl,dofmap,"comp_prof",&time); CHKERRA(ierr);

    if (ADVDIF_CHECK_JAC) {
      VOID_IT(argl);
      //  jobinfo_fields = "liq";
      //  argl.arg_add(AA_liq,PROFILE|PFMAT);
      argl.arg_add(AA_stage[kstage],PROFILE|PFMAT);
      ierr = assemble(mesh,argl,dofmap,"comp_prof",&time); CHKERRA(ierr);
    }
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  ierr = opt_read_vector(mesh,x,dofmap,MY_RANK); CHKERRA(ierr);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Hook stuff
  HookList hook_list;
  hook_list.init(*mesh,*dofmap,advdif_hook_factory);

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
    hook_list.time_step_pre(time_star.time(),tstep);

    for (int inwt=0; inwt<nnwt; inwt++) {
      
      
      for (int kstage=0; kstage<nstage; kstage++) {
	
	scal=0;
	ierr = VecSet(res_out,scal); CHKERRA(ierr);
	ierr = VecSet(dx_out,scal); CHKERRA(ierr);
	
	if (nstage==1) {
	  jobinfo_fields = "gasliq";
	  nnwt_in = nnwt_liq;
	  omega_newton_in = omega_newton_liq;
	}
	
	else if (kstage==0) {
	  jobinfo_fields = "liq";
	  nnwt_in = nnwt_liq;
	  omega_newton_in = omega_newton_liq;
	}
	else if (kstage==1) {
	  jobinfo_fields = "gas";
	  nnwt_in = nnwt_gas;
	  omega_newton_in = omega_newton_gas;
	}
	else if (kstage==2) {
	  jobinfo_fields = "kep";
	  nnwt_in = nnwt_liq;
	  omega_newton_in = omega_newton_liq;
	}
	
	
	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
	// Stages loop  "kstage"
	
	for (int inwt_in=0; inwt_in<nnwt_in; inwt_in++) {
	  
	  // Initializes res
	  scal=0;
	  ierr = VecSet(res,scal); CHKERRA(ierr);
	  ierr = VecSet(dx,scal); CHKERRA(ierr);
	  
	  if (comp_mat_each_time_step_g) {
	    
	    ierr = A_stage[kstage]->clean_mat(); CHKERRA(ierr);
            if (ADVDIF_CHECK_JAC) {
              ierr = AA_stage[kstage]->clean_mat(); CHKERRA(ierr);
            }
	    VOID_IT(argl);
	    argl.arg_add(&xold,IN_VECTOR);
            if (!ADVDIF_CHECK_JAC) argl.arg_add(&x,IN_VECTOR);
	    else argl.arg_add(&x,PERT_VECTOR);
	    argl.arg_add(&res,OUT_VECTOR);
	    argl.arg_add(&dtmin,VECTOR_MIN);
	    argl.arg_add(A_stage[kstage],OUT_MATRIX|PFMAT);
	    argl.arg_add(&glob_param,USER_DATA);
            if (ADVDIF_CHECK_JAC)
              argl.arg_add(AA_stage[kstage],OUT_MATRIX_FDJ|PFMAT);
	    
	    if (measure_performance) {
	      ierr = measure_performance_fun(mesh,argl,dofmap,"comp_res",
					     &time_star); CHKERRA(ierr);
	      PetscFinalize();
	      exit(0);
	    }
	    
	    // jobinfo_fields = "gas";
	    
	    ierr = assemble(mesh,argl,dofmap,"comp_res",&time_star); CHKERRA(ierr);
	    debug2.trace("After residual computation.");
	    
	    if (!print_linear_system_and_stop || solve_system) {
	      ierr = A_stage[kstage]->solve(res,dx); CHKERRA(ierr);
	    }
	    
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

	    // jobinfo_fields = "gas";
	    ierr = assemble(mesh,argl,dofmap,"comp_res",&time_star); CHKERRA(ierr);
	    
	    if (!print_linear_system_and_stop || solve_system) {
	      ierr =A_stage[kstage]->solve(res,dx); CHKERRA(ierr);
	    }
	  }
	  

	if (print_linear_system_and_stop &&
	    inwt==inwt_stop && tstep==time_step_stop) {
	  PetscPrintf(PETSCFEM_COMM_WORLD,
		      "Printing residual and matrix for debugging and stopping..\n");
	  PetscViewer matlab;
	  ierr = PetscViewerASCIIOpen(PETSCFEM_COMM_WORLD,
				      "mat.output",&matlab); CHKERRA(ierr);
	  ierr = PetscViewerSetFormat_WRAPPER(matlab,
					      PETSC_VIEWER_ASCII_MATLAB,"res");
	  ierr = VecView(res,matlab);
	  if (solve_system) {
	    ierr = PetscViewerSetFormat_WRAPPER(matlab,
						PETSC_VIEWER_ASCII_MATLAB,"dx");
	    ierr = VecView(dx,matlab);
	  }
	  ierr = PetscViewerSetFormat_WRAPPER(matlab,
					      PETSC_VIEWER_ASCII_MATLAB,"A");
	  ierr = A_stage[kstage]->view(matlab);
	  print_vector(save_file_res.c_str(),res,dofmap,&time); // debug:=
          if (ADVDIF_CHECK_JAC) {
            ierr = PetscViewerSetFormat_WRAPPER(matlab, 
                                                PETSC_VIEWER_ASCII_MATLAB,"AA");
            ierr = AA_stage[kstage]->view(matlab);
          }
	  PetscFinalize();
	  exit(0);
	}

	if (inwt_in==0){
	  scal=1;
	  ierr = VecAXPY(res_out,scal,res);
	  ierr = VecAXPY(dx_out,scal,dx);
	}

	ierr  = VecNorm(res,NORM_2,&normres); CHKERRA(ierr);
	PetscPrintf(PETSCFEM_COMM_WORLD,
		    "Newton subiter (inner) %d, stage  %d, norm_res  = %10.3e\n",
		    inwt_in,kstage,normres);
	scal=omega_newton_in/alpha;
	ierr = VecAXPY(x,scal,dx);
	if (normres < tol_newton) break;

	}  // end of inwt (inner) loop

	} // end of kstage loop

      ierr = VecCopy(res_out,res);

      ierr  = VecNorm(res,NORM_2,&normres); CHKERRA(ierr);
      if (inwt==0) normres_ext = normres;
      PetscPrintf(PETSCFEM_COMM_WORLD,
		  "Newton subiter (outer) %d, norm_res  = %10.3e\n",
		  inwt,normres);
      //      scal=omega_newton;
      if (normres < tol_newton) break;
    }  // end of inwt (outer) loop

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
    ierr = VecAXPY(dx,scal,xold);
    ierr  = VecNorm(dx,NORM_2,&delta_u); CHKERRA(ierr);

    if (tstep % nsave == 0){
      PetscPrintf(PETSCFEM_COMM_WORLD,
		  " --------------------------------------\n"
		  "Time step: %d\n"
		  " --------------------------------------\n",
		  tstep);

      print_vector(save_file.c_str(),x,dofmap,&time);
      if (print_residual)
	print_vector(save_file_res.c_str(),res,dofmap,&time);
    }

    if (ngather>0) {
      gather_values.resize(ngather,0.);
      for (int j=0; j<ngather; j++) gather_values[j] = 0.;
      arglf.clear();
      arglf.arg_add(&x,IN_VECTOR);
      arglf.arg_add(&xold,IN_VECTOR);
      arglf.arg_add(&gather_values,VECTOR_ADD);
      ierr = assemble(mesh,arglf,dofmap,"gather",&time);
      CHKERRA(ierr);
    }

    hook_list.time_step_post(time_star.time(),tstep,gather_values);

    if (ngather>0) {
      // Print gathered values
      if (MY_RANK==0) {
	if (gather_file == "") {
	  printf("Gather results: \n");
	  for (unsigned int j=0; j < gather_values.size(); j++) 
	    printf("v_component_%d = %12.10e\n",j,gather_values[j]);
	} else {
	  gather_file_f = fopen(gather_file.c_str(),"a");
	  for (unsigned int j=0; j<gather_values.size(); j++) 
	    fprintf(gather_file_f,"%12.10e ",gather_values[j]);
	  fprintf(gather_file_f,"\n");
	  fclose(gather_file_f);
	}
      }
    }

    PetscPrintf(PETSCFEM_COMM_WORLD,
		"time_step %d, time: %g, delta_u = %10.3e\n",
		tstep,time_,delta_u);
    print_vector_rota(save_file_pattern.c_str(),x,dofmap,
		      &time,tstep-1,nsaverot,nrec,nfile);

    if (print_some_file!="" && tstep % nsome == 0)
      print_some(save_file_some.c_str(),x,dofmap,node_list,&time);

    if (normres_ext < tol_steady) break;

  }
  hook_list.close();

  print_vector(save_file.c_str(),x,dofmap,&time);
  if (print_residual)
    print_vector(save_file_res.c_str(),res,dofmap,&time);
  if (report_option_access && MY_RANK==0) TextHashTable::print_stat();

  ierr = VecDestroy(&x); CHKERRA(ierr);
  ierr = VecDestroy(&xold); CHKERRA(ierr);
  ierr = VecDestroy(&dx); CHKERRA(ierr);
  ierr = VecDestroy(&res); CHKERRA(ierr);
  ierr = VecDestroy(&res_out); CHKERRA(ierr);
  ierr = VecDestroy(&dx_out); CHKERRA(ierr);

    for (int kstage=0; kstage<nstage; kstage++) {

#ifdef DIAG_MAT_MATRIX
//  ierr = MatDestroy(&A_liq); CHKERRA(ierr);
//  ierr = MatDestroy(&A_gas); CHKERRA(ierr);
  ierr = MatDestroy(&A_stage[kstage]); CHKERRA(ierr);
#endif

//  delete A_liq,A_gas;
//  delete AA_liq,AA_gas;
  delete A_stage[kstage];
  delete AA_stage[kstage];
}

  PetscFinalize();
  exit(0);
}
