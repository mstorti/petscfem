//__INSERT_LICENSE__
//$Id: advdif.cpp,v 1.71.10.2 2007/02/23 19:18:07 rodrigop Exp $

#include <src/debug.h>
#include <set>

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/pfmat.h>
#include <src/hook.h>

#include <src/fm2stats.h>
#include <src/fastlib2.h>
#include <src/fm2prod.h>

#include "advective.h"

#include <time.h>

static char help[] = "Basic finite element program.\n\n";
int abso_inwt;

extern int print_internal_loop_conv_g;
extern int consistent_supg_matrix_g;
extern int  local_time_step_g;
extern int  comp_mat_each_time_step_g;

#define VECVIEW(name,label) \
ierr = PetscViewerSetFormat(matlab, \
		       PETSC_VIEWER_ASCII_MATLAB,#label); \
ierr = VecView(name,matlab); CHKERRA(ierr)

// PETSc now doesn't have the string argument that represents the variable name
// so that I will use this wrapper until I find how to set names in Ascii matlab viewers.
#define PetscViewerSetFormat_WRAPPER(viewer,format,name) \
          PetscViewerSetFormat(viewer,format)

int bubbly_main();
int fsi_main();
int bdf_main();
int dual_time_main();

Hook *advdif_hook_factory(const char *name);

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int VecIsFinite(Vec x) {
  double *xp;
  int ierr,size;
  ierr = VecGetArray(x,&xp); CHKERRQ(ierr);
  ierr = VecGetLocalSize(x,&size); CHKERRQ(ierr);
  int locok = 1, ok;
  for (int j=0; j<size; j++) 
    if (!isfinite(xp[j])) { locok = 0; break; }
  ierr = VecRestoreArray(x,&xp);
  ierr = MPI_Allreduce(&locok,&ok,1,MPI_INT,
                       MPI_LAND,PETSCFEM_COMM_WORLD);
  return ok;
}
#endif

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
#undef __FUNC__
#define __FUNC__ "advdif_main"
int advdif_main(int argc,char **args) {

#ifdef DO_SIZE_STATS
  prod2_subcache_t::do_size_stats=1;
#endif
  // prod2_subcache_t::nmax=0;
  PetscFemInitialize(&argc,&args,(char *)0,help);
  
#define CNLEN 100
  PetscBool flg;
  char code_name[CNLEN];
  int ierr = PetscOptionsGetString(PETSC_NULL,"-code",code_name,CNLEN,&flg);

  if (flg) {
    if (!strcmp(code_name,"fsi")) return fsi_main();
    if (!strcmp(code_name,"bdf")) return bdf_main();
    if (!strcmp(code_name,"bubbly")) return bubbly_main();
    if (!strcmp(code_name,"dual_time")) return dual_time_main();
    PETSCFEM_ERROR("Unknown -code option: \"%s\"\n",code_name);
  }
  
  Vec     x, dx, xold, res; /* approx solution, RHS, residual*/
  PFMat *A,*AA;			// linear system matrix 
  PFMat *A_tet, *A_tet_c;
  double  *sol, scal, normres, normres_ext=NAN;    /* norm of solution error */
  int     i, n = 10, col[3], its, size, node,
    jdof, k, kk, nfixa,
    kdof, ldof, lloc, ndim, nel, nen, neq, nu,
    myrank;
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
	      "-------- Generic Advective-Diffusive  module ---------\n");

  Debug debug(0,PETSCFEM_COMM_WORLD);
  GLOBAL_DEBUG = &debug;

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

  //o Activate debugging
  GETOPTDEF(int,activate_debug,0);
  if (activate_debug) {
    debug.activate();
    Debug::init();
  }
  //o Activate printing in debugging
  GETOPTDEF(int,activate_debug_print,0);
  if (activate_debug_print) debug.activate("print");
  //o Activate report of memory usage
  GETOPTDEF(int,activate_debug_memory_usage,0);
  if (activate_debug_memory_usage) debug.activate("memory_usage");

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
  //o Measure performance of the  #comp_mat_res#  jobinfo. 
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
  //o Solve system before  #print\_linear_system_and_stop# 
  GETOPTDEF(int,solve_system,1);
  //o If #print_linear_system_and_stop# is active,
  // then print system in this Newton iteration 
  GETOPTDEF(int,inwt_stop,0);
  //o If #print_linear_system_and_stop# is active,
  // then print system in this time step
  GETOPTDEF(int,time_step_stop,1);
  //o Print the residual each  #nsave#  steps. 
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
  //o Update jacobian each $n$-th time step. 
  GETOPTDEF(int,update_jacobian_steps,0);
  //o Use IISD (Interface Iterative Subdomain Direct) or not.
  GETOPTDEF(int,use_iisd,0);
  //o Type of solver. May be  #iisd#  or  #petsc# . 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,solver,petsc);
  if (use_iisd) solver = string("iisd");

  // Use IISD (Interface Iterative Subdomain Direct) or not.
  // A_tet = (use_iisd ? &IISD_A_tet : &PETSc_A_tet);
  A  = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());
  AA = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());

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
  //o Flag if steady solution or not (uses Dt=inf). If  #steady# 
  // is set to 1, then the computations are as if $\Dt=\infty$. 
  // The value of  #Dt#  is used for printing etc... If  #Dt# 
  // is not set and  #steady#  is set then  #Dt#  is set to one.
  GETOPTDEF(int,steady,0);
  if (steady && Dt==0.) Dt=1.;
  glob_param.Dt = Dt;
  glob_param.steady = steady;
  //o The parameter of the trapezoidal rule
  // for temporal integration. 
  GETOPTDEF(double,alpha,1.);
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

  //o Flag for launching RENORM process
  GETOPTDEF(int,RENORM_flag,0);

  //o Flag for volume conservation in ALE
  GETOPTDEF(int,VCONS_ALE,0);
  //o Total initial volume for reference - volume conservation - ALE FS
  TGETOPTDEF(GLOBAL_OPTIONS,double,initial_volume,0.);
  //o Free surface horizontal projection - volume conservation - ALE FS
  TGETOPTDEF(GLOBAL_OPTIONS,double,efs,0.);
  //o Volume control coefficient - volume conservation - ALE FS
  TGETOPTDEF(GLOBAL_OPTIONS,double,Cvol_ale,0.1);

  comp_mat_each_time_step_g = 
    consistent_supg_matrix || local_time_step;

  //o Counts time from here.
  GETOPTDEF(double,start_comp_time,0.);
  //o Counts time from here. (superseded by
  //  #start_comp_time# for compatibility with Navier-Stokes
  //  module). 
  GETOPTDEF(double,start_time,NAN);
  if (!isnan(start_time) && start_comp_time!=0.0)
    start_comp_time = start_time;
  
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
  //o Computes jacobian of residuals and prints to a file.
  //  May serve to debug computation of the analytic jacobians. 
  TGETOPTDEF_ND(mesh->global_options,int,verify_jacobian_with_numerical_one,0);
  //o Check whether the states are finite (not Inf or NaN).
  //  If this happens the program stops. 
  GETOPTDEF(int,check_for_inf,0);

#define INF INT_MAX
  //o Update jacobian each $n$-th time step. 
  GETOPTDEF(int,update_jacobian_start_steps,INF);
  //o Update jacobian only until n-th Newton subiteration. 
  // Don't update if null. 
  GETOPTDEF(int,update_jacobian_iters,1);
  PETSCFEM_ASSERT0(update_jacobian_iters>=1,"Out of range");  
  //o Update jacobian each $n$-th Newton iteration
  GETOPTDEF(int,update_jacobian_start_iters,INF);
  PETSCFEM_ASSERT0(update_jacobian_start_iters>=0,"Out of range");  
#undef INF

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
  time.set(start_comp_time);
  glob_param.time = &time;

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

  // Set pointers in glob_param
  glob_param.x = x;
  glob_param.xold = xold;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // initialize state vectors
  scal=0;
  ierr = VecSet(x,scal); CHKERRA(ierr);

  arg_list argl,arglf;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Hook stuff
  HookList hook_list;
  hook_list.init(*mesh,*dofmap,advdif_hook_factory);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Compute  profiles
  debug.trace("Computing profiles...");
  VOID_IT(argl);
  argl.arg_add(A,PROFILE|PFMAT);
  ierr = assemble(mesh,argl,dofmap,"comp_prof",&time); CHKERRA(ierr);
  debug.trace("After computing profile.");

  if (ADVDIF_CHECK_JAC) {
    VOID_IT(argl);
    argl.arg_add(AA,PROFILE|PFMAT);
    ierr = assemble(mesh,argl,dofmap,"comp_prof",&time); CHKERRA(ierr);
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  ierr = opt_read_vector(mesh,x,dofmap,MY_RANK); CHKERRA(ierr);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // This is for taking statistics of the
  // CPU time consumed by a time steptime
  Chrono chrono; 
#define STAT_STEPS 5
  double cpu[STAT_STEPS],cpuav;
  int update_jacobian_this_step,update_jacobian_this_iter;
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

    if (RENORM_flag)
      ierr = read_vector("state-ren.tmp",x,dofmap,MY_RANK); CHKERRA(ierr);

    // Volume preservation in ALE
    if ((VCONS_ALE) && (tstep>2)){
      FILE *fd; // File pointer
      char buff[1024];
      
      if ((fd = fopen("gather.out", "r")) != NULL) // open file
	{
	  fseek(fd, -18, SEEK_END); 
	  //	  fseek(fd, 0, SEEK_SET); // make sure start from 0
	  
	  while(!feof(fd))
	    {
	      memset(buff, 0x00, 1024); // clean buffer
	      fscanf(fd, "%[^\n]\n", buff); // read file 
	    }
	  printf("Current volume :: %s\n", buff);
	  double current_volume = atof(buff);
	  double delta_eta_V = (current_volume - initial_volume)/efs;
	  ierr = VecShift(x,-Cvol_ale*delta_eta_V); CHKERRA(ierr);
	  PetscPrintf(PETSCFEM_COMM_WORLD,"Delta eta %.10f \n",-0.1*delta_eta_V);
	}
      fclose(fd);
    }
    ierr = VecCopy(x,xold);
    //    hook_list.time_step_pre(time_star.time(),tstep);
    hook_list.time_step_pre(time.time()+Dt,tstep); //hook needs t_{n+1}

    for (int inwt=0; inwt<nnwt; inwt++) {
      abso_inwt = inwt;

      // Initializes res
      scal=0;
      ierr = VecSet(res,scal); CHKERRA(ierr);

      if (comp_mat_each_time_step_g) {

	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
	// ierr = A->build_ksp(GLOBAL_OPTIONS); CHKERRA(ierr); 

	ierr = A->clean_mat(); CHKERRA(ierr); 
        if (ADVDIF_CHECK_JAC) {
          ierr = AA->clean_mat(); CHKERRA(ierr);
        }
	VOID_IT(argl);
	argl.arg_add(&xold,IN_VECTOR);

        if (!ADVDIF_CHECK_JAC) argl.arg_add(&x,IN_VECTOR);
        else argl.arg_add(&x,PERT_VECTOR);

	argl.arg_add(&res,OUT_VECTOR);
	argl.arg_add(&dtmin,VECTOR_MIN);
	argl.arg_add(A,OUT_MATRIX|PFMAT);
	argl.arg_add(&glob_param,USER_DATA);
        if (ADVDIF_CHECK_JAC)
          argl.arg_add(AA,OUT_MATRIX_FDJ|PFMAT);

	if (measure_performance) {
	  ierr = measure_performance_fun(mesh,argl,dofmap,"comp_res",
					 &time_star); CHKERRA(ierr);
	  PetscFinalize();
	  exit(0);
	}
	debug.trace("Before residual computation...");
	ierr = assemble(mesh,argl,dofmap,"comp_res",&time_star); CHKERRA(ierr);
	debug.trace("After residual computation.");

	if (!print_linear_system_and_stop || solve_system) {
	  debug.trace("Before solving linear system...");
	  ierr = A->solve(res,dx); CHKERRA(ierr); 
	  debug.trace("After solving linear system.");
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
	ierr = assemble(mesh,argl,dofmap,"comp_res",&time_star); CHKERRA(ierr);

	if (!print_linear_system_and_stop || solve_system) {
	  ierr = A->solve(res,dx); CHKERRA(ierr); 
	}
      }

      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
      // FEM matrix jacobian debug  (perturbation)
#if 0
	PetscViewer matlab;
	if (verify_jacobian_with_numerical_one) {
	  ierr = PetscViewerASCIIOpen(PETSCFEM_COMM_WORLD,
				      "system.dat.tmp",&matlab); CHKERRA(ierr);
	  ierr = PetscViewerSetFormat_WRAPPER(matlab,
					      PETSC_VIEWER_ASCII_MATLAB,
					      "atet"); CHKERRA(ierr);
	  
	  ierr = A_tet->view(matlab); CHKERRQ(ierr); 
	  
	  ierr = A_tet_c->duplicate(MAT_DO_NOT_COPY_VALUES,*A_tet); CHKERRA(ierr);
	  ierr = A_tet->clean_mat(); CHKERRA(ierr); 
	  ierr = A_tet_c->clean_mat(); CHKERRA(ierr); 
	  
	  argl.clear();
	  argl.arg_add(&x,PERT_VECTOR);
	  argl.arg_add(&xold,IN_VECTOR);
	  argl.arg_add(A_tet_c,OUT_MATRIX_FDJ|PFMAT);
	  
	  argl.arg_add(A_tet,OUT_MATRIX|PFMAT);
	  argl.arg_add(&hmin,VECTOR_MIN);

	  argl.arg_add(&glob_param,USER_DATA);
	  argl.arg_add(wall_data,USER_DATA);
	  ierr = assemble(mesh,argl,dofmap,jobinfo,
			  &time_star); CHKERRA(ierr);
	  
	  ierr = PetscViewerSetFormat_WRAPPER(matlab,
					      PETSC_VIEWER_ASCII_MATLAB,"atet_fdj"); CHKERRA(ierr);
	  ierr = A_tet_c->view(matlab); CHKERRQ(ierr); 
	  
	  PetscFinalize();
	  exit(0);
	}
#endif
      
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
	ierr = A->view(matlab);
	print_vector(save_file_res.c_str(),res,dofmap,&time); // debug:=
        if (ADVDIF_CHECK_JAC) {
          ierr = PetscViewerSetFormat_WRAPPER(matlab, 
                                              PETSC_VIEWER_ASCII_MATLAB,"AA");
          ierr = AA->view(matlab);
        }
	PetscFinalize();
	exit(0);
      }

      ierr  = VecNorm(res,NORM_2,&normres); CHKERRA(ierr);
      if (check_for_inf) {
#if 0
        PETSCFEM_ASSERT0(VecIsFinite(res),
                         "Detected Inf or NaN values in residual vector");  
#else
        PETSCFEM_ASSERT0(isfinite(normres),
                         "Detected Inf or NaN values in residual vector");  
#endif
      }
      if (inwt==0) normres_ext = normres;
      PetscPrintf(PETSCFEM_COMM_WORLD,
		  "Newton subiter %d, norm_res  = %10.3e\n",
		  inwt,normres);
      scal=omega_newton/alpha;
      ierr = VecAXPY(x,scal,dx);
      if (check_for_inf) {
#if 0
        PETSCFEM_ASSERT0(VecIsFinite(x),
                         "Detected Inf or NaN values in state vector");  
#else
        double normx;
        ierr  = VecNorm(x,NORM_2,&normx); CHKERRA(ierr);
        PETSCFEM_ASSERT0(isfinite(normx),
                         "Detected Inf or NaN values in state vector");  
#endif
      }
      if (normres < tol_newton) break;
      } // end of Newton loop

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

    //    hook_list.time_step_post(time_star.time(),tstep,gather_values);
    hook_list.time_step_post(time.time()+Dt,tstep,gather_values);

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
#ifdef DIAG_MAT_MATRIX
  ierr = MatDestroy(&A); CHKERRA(ierr); 
#endif
  
  delete A;
  delete AA;

#ifdef DO_SIZE_STATS
  prod2_subcache_t::report_stats();
#endif

  PetscFinalize();
  exit(0);
}
