//__INSERT_LICENSE__
//$Id$
#include <src/debug.h>
#include <malloc.h>

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/sttfilter.h>
#include <src/pfmat.h>
#include <src/hook.h>

// PETSc now doesn't have the string argument that represents the variable name
// so that I will use this wrapper until I find how to set names in Ascii matlab viewers.
#define PetscViewerSetFormat_WRAPPER(viewer,format,name) \
          PetscViewerSetFormat(viewer,format)

#include <applications/ns/nsi_tet.h>
static char help[] = "PETSc-FEM Navier Stokes module\n\n";

extern int MY_RANK,SIZE;
extern double mmv_delta,mmv_d2fd,mmv_dfd;
extern double min_quality,min_volume;
extern double mmv_functional;
extern int    tarea,tangled_mesh;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Creates hooks depending on the name. 
    @param name (input) the name of the hook. 
    @return a pointer to the hook. */ 
extern Hook *ns_hook_factory(const char *name);

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
#undef __FUNC__
extern GlobParam *GLOB_PARAM;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define __FUNC__ "struct_main"
int mmove2_main() {

  // approx solution, RHS, residual
  Vec x, dx, xold, dx_step, res, res_delta, resp; 
  PetscViewer matlab;
  PFMat *Ap=NULL;	// linear system matrix 
  double  norm, *sol, scal;	// norm of solution error
  double  mmv_delta_increment,resdelta_dot_dx;
  double  min_quality_old;
  double  local_mmv_dfd,local_mmv_d2fd;
  double  mmv_functional_g,mmv_dfd_g,mmv_d2fd_g;
  double  min_quality_g,min_volume_g;
  int     ierr, i, n = 10, size, node,
    jdof, k, kk, nfixa,
    kdof, ldof, lloc, nel, nen, neq, nu,
    myrank;
  PetscBool flg;
  // Initialize time
  Time time,time_old; 
  GlobParam glob_param;
  GLOB_PARAM = &glob_param;
  string save_file_res;
  BasicObject_application_factory = &BasicObject_ns_factory;
  
  // ierr = MatCreateShell(PETSC_COMM_WORLD,int m,int n,int M,int N,void *ctx,Mat *Ap)
  char fcase[FLEN+1],output_file[FLEN+1];
  Dofmap *dofmap;
  Mesh *mesh;
  // arglf:= argument list for computing gathered quantities as forces
  arg_list argl, arglf;
  vector<double> hmin;
  hmin.resize(1);

  print_copyright();
  PetscPrintf(PETSC_COMM_WORLD,
              "=============================================\n"
              "-------- MESH-MOVE VERSION 2 MODULE ---------\n"
              "=============================================\n"
              );

  Debug debug(0,PETSC_COMM_WORLD);
  GLOBAL_DEBUG = &debug;

  int activate_debug=0;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-activate_debug",
			    &activate_debug,&flg);

  if (activate_debug) {
    debug.activate();
    Debug::init();
  }
  ierr = PetscOptionsGetString(PETSC_NULL,"-case",fcase,FLEN,&flg);
  CHKERRA(ierr);
  if (!flg) {
    PetscPrintf(PETSC_COMM_WORLD,
		"Option \"-case <filename>\" not passed to PETSc-FEM!!\n");
    PetscFinalize();
    exit(0);
  }

  ierr = PetscOptionsGetString(PETSC_NULL,"-o",output_file,FLEN,&flg);
  CHKERRA(ierr);
  if (flg) { 
    PetscPrintf(PETSC_COMM_WORLD,"PETSc-FEM: NS module: "
		"redirecting output to \"%s\"\n",output_file);
    FILE *new_stdout = fopen(output_file,"w");
    if (!new_stdout) {
      PetscPrintf(PETSC_COMM_WORLD,"error redirecting output. "
		  "Couldn't open \"%s\"\n",output_file);
    } else {
      fclose(stdout);
      stdout = new_stdout;
    }
  }

  // Read the mesh
  debug.trace("Before readmesh.");
  ierr = read_mesh(mesh,fcase,dofmap,neq,SIZE,MY_RANK); CHKERRA(ierr); 
  debug.trace("After readmesh.");

  GLOBAL_OPTIONS = mesh->global_options;

  //o Time to start computations
  TGETOPTDEF(GLOBAL_OPTIONS,double,start_comp_time,0.);
  time.set(start_comp_time);
  // time_old.set(start_comp_time-Dt); // we should do this!!
  State state(x,time),state_old(xold,time_old);
  vector<State *> add_states;

  //o Additional states to be used by modules
  TGETOPTDEF(GLOBAL_OPTIONS,double,time_fac_epsilon,1e-3);
  assert(time_fac_epsilon>0.0);

#if 0
  //o If set, redirect output to this file.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,stdout_file,none);
  if (!strcmp(stdout_file.c_str(),"")) {
    fclose(stdout);
    stdout = fopen(stdout_file.c_str(),"w");
  }
#endif

  //o Additional states to be used by modules
  TGETOPTDEF(GLOBAL_OPTIONS,int,additional_states,0);
  //o Use fractional step or TET algorithm
  TGETOPTDEF_ND(GLOBAL_OPTIONS,int,fractional_step,0);
  //o Use fractional step or TET algorithm
  TGETOPTDEF_ND(GLOBAL_OPTIONS,int,reuse_mat,0);
  //o Fractional step uses symmetric matrices (only CG iterative KSP).
  TGETOPTDEF(GLOBAL_OPTIONS,int,fractional_step_use_petsc_symm,1);
  //o Solver combination for the fractional step method. May be #iisd#, 
  //  #lu#, #global_gmres#. 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,fractional_step_solver_combo,iisd);

  //o After computing the linear system for the
  //  predictor/momentum step print right hand side
  //  and solution vector, and stops.
  GETOPTDEF(int,stop_mom,0);
  //o After computing the linear system for the poisson step print
  // right hand side and solution vector, and stops. 
  GETOPTDEF(int,stop_poi,0);
  //o After computing the linear system for the projection step print
  // right hand side and solution vector, and stops. 
  GETOPTDEF(int,stop_prj,0);
  //o After computing the linear system for the
  //  predictor/momentum step print right hand side
  //  and solution vector, and stops.
  GETOPTDEF(int,stop_on_step,1);
  int do_stop;

  //o Activate debugging
  TGETOPTDEF_ND(GLOBAL_OPTIONS,int,activate_debug,0);
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

  HookList hook_list;
  hook_list.init(*mesh,*dofmap,ns_hook_factory);

  //o Dimension of the problem.
  GETOPTDEF(int,ndim,3);

  //o Scales displacement for ALE-like mesh relocation. 
  GETOPTDEF(double,displ_factor,0.1);
  //o Number of inner iterations for the global non-linear
  // Newton  problem. 
  GETOPTDEF(int,nnwt,1);
  //o Tolerance to solve the non-linear system (global Newton).
  GETOPTDEF(double,tol_newton,1e-8);
  //o 
  GETOPTDEF(int,use_mmv_predictor,0);
  //o 
  GETOPTDEF(double,mmv_h_star,0.75);
  //o Tolerance relative for the minimum quality in the Newton loop
  GETOPTDEF(double,tol_rel_quality,1.e-2);
  //o 
  GETOPTDEF(double,delta_fraction,0.);
  //o 
  GETOPTDEF(int,use_line_search_mmv,0);
  //o 
  GETOPTDEF(double,ls_relfac_fraction,2.);
  //o 
  GETOPTDEF(double,ls_eta_parameter,1.e-4);

  //o _T: vector<double>
  // _N: newton_relaxation_factor
  // _D: (none)
  // _DOC:
  //i_tex ../../doc/nsdoc.tex newton_relaxation_factor
  // _END
  vector<double> newton_relaxation_factor;
  mesh->global_options->get_entry("newton_relaxation_factor",
				  newton_relaxation_factor);
  // If entered empty vector, then set to 1.
  if (newton_relaxation_factor.size()==0) 
    newton_relaxation_factor.push_back(1.);
  // Check number of elements is odd
  assert(newton_relaxation_factor.size() % 2 ==1);
  // Check even entries are integer
  for (unsigned int j=1; j<newton_relaxation_factor.size(); j += 2) {
    double v = double(int(newton_relaxation_factor[j]));
    assert(v == newton_relaxation_factor[j]);
  }

  //o Computes jacobian of residuals and prints to a file.
  //  May serve to debug computation of the analytic jacobians. 
  GETOPTDEF(int,verify_jacobian_with_numerical_one,0);
  //o After computing the linear system solves it and prints Jacobian,
  // right hand side and solution vector, and stops. 
  GETOPTDEF(int,print_linear_system_and_stop,0);
  //o Print the residual each  #nsave#  steps. 
  GETOPTDEF(int,print_residual,0);
  //o Solve system before  #print_linear_system_and_stop# 
  GETOPTDEF(int,solve_system,1);
  //o Measure performance of the 'comp\_mat\_res' jobinfo. 
  GETOPTDEF(int,measure_performance,0);
  //o Check whether the states are finite (not Inf or NaN).
  //  If this happens the program stops. 
  GETOPTDEF(int,check_for_inf,0);

  //o Sets the save frequency in iterations 
  GETOPTDEF(int,nsave,10);
  //o Sets the frequency save for the ``rotary save'' mechanism. 
  //i_tex ../../doc/nsdoc.tex rotary_save
  GETOPTDEF(int,nsaverot,100);
  //o Sets the number of states saved in a given file
  // in the ``rotary save'' mechanism (see \ref{sec:rotary_save}
  GETOPTDEF(int,nrec,1);
  //o Sets the number of files in the ``rotary save'' mechanism. 
  // (see \ref{sec:rotary_save})
  GETOPTDEF(int,nfile,-1);

  //o Sets the save frequency in iterations for the ``print some''
  // mechanism. 
  //i_tex ../../doc/nsdoc.tex print_some
  GETOPTDEF(int,nsome,10000);
  //o The number of time steps. 
  GETOPTDEF(int,nstep,10000);
  //o The time step.
  GETOPTDEF(double,Dt,0.);
  //o Flag if steady solution or not (uses Dt=inf). If  #steady# 
  // is set to 1, then the computations are as if $\Dt=\infty$. 
  // The value of  #Dt#  is used for printing etc... If  #Dt# 
  // is not set and  #steady#  is set then  #Dt#  is set to one.
  GETOPTDEF(int,steady,0);
  if (steady && Dt==0.) Dt=1.;
  // Set values to be passed through global options
  glob_param.steady=steady;
  glob_param.Dt=Dt;
  glob_param.state = &state;
  glob_param.state_old = &state_old;
  //o Trapezoidal method parameter.  #alpha=1# :
  // Backward Euler.  #alpha=0# : Forward Euler.
  //  #alpha=0.5# : Crank-Nicholson. 
  GETOPTDEF(double,alpha,1.); 
  glob_param.alpha=alpha;

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

  //o Use IISD (Interface Iterative Subdomain Direct) or not.
  GETOPTDEF(int,use_iisd,0);
  //o Type of solver. May be  #iisd#  or  #petsc# . 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,solver,petsc);
  if (use_iisd) solver = string("iisd");
  //o Type of solver for the projection and momentum steps
  // (fractional-step). May be  #iisd#  or  #petsc# .
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,solver_mom,petsc);
  if (use_iisd) solver_mom = string("iisd");

  //o The pattern to generate the file name to save in for
  // the rotary save mechanism.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file_pattern,outvector%d.out);

  //o The name of the file to save the state vector. 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file,outvector.out);
  save_file_res = save_file + string(".res");

  //o Name of file where to read the nodes for the ``print some'' 
  // feature. 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,print_some_file,<none>);
  //o Name of file where to save node values for the ``print some'' 
  // feature. 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file_some,outvsome.out);
  //o Access mode to the ``some'' file. If 0 rewind file. If 1 
  //  append to previous  results.
  TGETOPTDEF(GLOBAL_OPTIONS,int,save_file_some_append,1);
  //o Print, after execution, a report of the times a given option
  // was accessed. Useful for detecting if an option was used or not.
  GETOPTDEF(int,report_option_access,1);

  if (print_some_file=="<none>")
    print_some_file = "";
  set<int> node_list;
  print_some_file_init(mesh->global_options,
		       print_some_file.c_str(),
		       save_file_some.c_str(),node_list,
		       save_file_some_append);

  // initialize vectors
  dofmap->create_MPI_vector(x);
  State sx(x,time);		// convert to State

  Ap = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());

  ierr = VecDuplicate(x,&xold); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx_step); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx); CHKERRA(ierr);
  ierr = VecDuplicate(x,&res); CHKERRA(ierr);
  ierr = VecDuplicate(x,&res_delta); CHKERRA(ierr);
  ierr = VecDuplicate(x,&resp); CHKERRA(ierr);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // initialize state vectors
  scal=0;
  ierr = VecSet(x,scal); CHKERRA(ierr);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Compute  profiles
  debug.trace("Computing profiles...");
  argl.clear();
  argl.arg_add(Ap,PROFILE|PFMAT,"A");
  ierr = assemble(mesh,argl,dofmap,"comp_prof",&time); CHKERRA(ierr); 
  debug.trace("After computing profile.");

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  debug.trace("Before reading initial vector...");
  ierr = opt_read_vector(mesh,x,dofmap,MY_RANK); CHKERRA(ierr);
  debug.trace("After reading initial vector...");

  // update_jacobian_this_step:= Flags whether this step the
  // jacobian should be updated or not 
  int update_jacobian_this_step,
    update_jacobian_this_iter, tstep_start=1;
  for (int tstep=tstep_start; tstep<=nstep; tstep++) {

    hook_list.time_step_pre(time.time(),tstep);
    
    //#define MMV_DBG
#ifdef MMV_DBG
    printf("x prev to project\n");
    ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);
    CHKERRA(ierr);
#endif

    // Inicializacion del paso
    ierr = VecCopy(x,xold);
    ierr = VecCopy(xold,dx_step);

    time_old.set(time.time());

    if (use_mmv_predictor){
      // Computes a better starting point based on the solution
      // of a differential problem
      // res = RES(u^n,t^n);
      tarea = 1; // Tarea = 1 calcula residuo y matriz
      mmv_functional = mmv_dfd = mmv_d2fd = 0.;
      min_quality = 1.;
      min_volume  = 1.e10;
      argl.clear();
      state.set_time(time);
      state_old.set_time(time_old);
      argl.arg_add(&state,IN_VECTOR|USE_TIME_DATA,"state");
      argl.arg_add(&state_old,IN_VECTOR|USE_TIME_DATA,"state_old");
      argl.arg_add(&res,OUT_VECTOR,"res");
      argl.arg_add(&res_delta,OUT_VECTOR,"res_delta");
      argl.arg_add(Ap,OUT_MATRIX|PFMAT,"A");
      argl.arg_add(&hmin,VECTOR_MIN,"hmin");
      argl.arg_add(&glob_param,USER_DATA,"glob_param");

      scal=0;
      ierr = VecSet(res,scal); CHKERRA(ierr);
      ierr = Ap->clean_mat(); CHKERRA(ierr); 
      debug.trace("Before residual computation...");
      ierr = assemble(mesh,argl,dofmap,"comp_mat_res",&time);
      CHKERRA(ierr);
      debug.trace("After residual computation.");

      MPI_Allreduce(&mmv_functional, &mmv_functional_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&mmv_dfd, &mmv_dfd_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&mmv_d2fd, &mmv_d2fd_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&min_quality, &min_quality_g, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&min_volume, &min_volume_g, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      
      // res = RES(u^n,t^n+epsilon);
      mmv_functional = mmv_dfd = mmv_d2fd = 0.;
      min_quality = 1.;
      min_volume  = 1.e10;
      double epsilon = time_fac_epsilon*Dt;
      time.inc(epsilon);
      argl.clear();
      state.set_time(time);
      state_old.set_time(time_old);
      argl.arg_add(&state,IN_VECTOR|USE_TIME_DATA,"state");
      argl.arg_add(&state_old,IN_VECTOR|USE_TIME_DATA,"state_old");
      argl.arg_add(&resp,OUT_VECTOR,"res");
      argl.arg_add(&res_delta,OUT_VECTOR,"res_delta");
      argl.arg_add(Ap,OUT_MATRIX|PFMAT,"A");
      argl.arg_add(&hmin,VECTOR_MIN,"hmin");
      argl.arg_add(&glob_param,USER_DATA,"glob_param");
      
      scal=0;
      ierr = VecSet(resp,scal); CHKERRA(ierr);
      ierr = Ap->clean_mat(); CHKERRA(ierr); 
      debug.trace("Before residual computation...");
      ierr = assemble(mesh,argl,dofmap,"comp_mat_res",&time);
      CHKERRA(ierr);
      debug.trace("After residual computation.");

      MPI_Allreduce(&mmv_functional, &mmv_functional_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&mmv_dfd, &mmv_dfd_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&mmv_d2fd, &mmv_d2fd_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&min_quality, &min_quality_g, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&min_volume, &min_volume_g, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      
#ifdef MMV_DBG
      printf("resp: ");
      ierr = VecView(resp,PETSC_VIEWER_STDOUT_WORLD);
      CHKERRA(ierr);
#endif

      // resp = (resp-res)
      scal = -1.;
      ierr = VecAXPY(resp,scal,res); CHKERRQ(ierr); 
    
#ifdef MMV_DBG
      printf("(resp-res)/time_fac_epsilon: ");
      ierr = VecView(resp,PETSC_VIEWER_STDOUT_WORLD);
      CHKERRA(ierr);
#endif

      scal = 1./time_fac_epsilon;
      ierr = VecScale(resp,scal);
      ierr = Ap->solve(resp,dx); CHKERRA(ierr); 
      
#ifdef MMV_DBG
      ierr = Ap->view(matlab); CHKERRA(ierr); 
      printf("dx: ");
      ierr = VecView(dx,PETSC_VIEWER_STDOUT_WORLD);
      CHKERRA(ierr);
#endif

      // x = x+dx
      scal = 1.0;
      ierr = VecAXPY(x,scal,dx); CHKERRA(ierr); 
    }

    time.set(time_old.time());
    time.inc(Dt);
    // print_vector("state-proj.tmp",x,dofmap,&time);

    // Esto hay que mejorarlo, es solo una propuesta 
    // para inicializar el delta
    tarea = 0; // Tarea = 0 solo calcula la calidad, 
               // no calcula residuo y matriz
    mmv_functional = mmv_dfd = mmv_d2fd = 0.;
    min_quality = 1.;
    min_volume  = 1.e10;
    argl.clear();
    state.set_time(time);
    state_old.set_time(time_old);
    argl.arg_add(&state,IN_VECTOR|USE_TIME_DATA,"state");
    argl.arg_add(&state_old,IN_VECTOR|USE_TIME_DATA,"state_old");
    argl.arg_add(&resp,OUT_VECTOR,"res");
    argl.arg_add(&res_delta,OUT_VECTOR,"res_delta");
    argl.arg_add(Ap,OUT_MATRIX|PFMAT,"A");
    argl.arg_add(&hmin,VECTOR_MIN,"hmin");
    argl.arg_add(&glob_param,USER_DATA,"glob_param");

    scal=0;
    ierr = VecSet(resp,scal); CHKERRA(ierr);
    ierr = Ap->clean_mat(); CHKERRA(ierr); 
    debug.trace("Before residual computation...");
    ierr = assemble(mesh,argl,dofmap,"comp_mat_res",&time);
    CHKERRA(ierr);
    debug.trace("After residual computation.");

    MPI_Allreduce(&mmv_functional, &mmv_functional_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&mmv_dfd, &mmv_dfd_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&mmv_d2fd, &mmv_d2fd_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&min_quality, &min_quality_g, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&min_volume, &min_volume_g, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    min_quality_old = min_quality_g;

    // tangled_mesh = (min_volume > 0. ? 0 : 1);
    mmv_delta    = (min_volume_g > 1.e-15 ? 0. : \
		    0.25*mmv_h_star*min_volume_g/(pow(mmv_h_star,2.)-1.)+1.e-10);
    tangled_mesh = (mmv_delta <= 1.e-15 ? 0 : 1);
    PetscPrintf(PETSC_COMM_WORLD,"delta: %12.10f\n",mmv_delta);
    PetscPrintf(PETSC_COMM_WORLD,"min vol: %12.10f\n",min_volume_g);

    if (!MY_RANK) printf("Time step: %d, time: %g %s\n",
			 tstep,time.time(),(steady ? " (steady) " : ""));

    //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    // NEWTON-RAPHSON ALGORITHM
    
    double normres_external=NAN;
    for (int inwt=0; inwt<nnwt; inwt++) {
      
      glob_param.inwt = inwt;
      // Initialize step
      
      tarea = 1;
      mmv_functional = mmv_dfd = mmv_d2fd = 0.;
      min_quality = 1.;
      min_volume  = 1.e10;
      
      ierr = Ap->clean_mat(); CHKERRA(ierr); 

      scal=0;
      ierr = VecSet(res,scal); CHKERRA(ierr);
      ierr = VecSet(res_delta,scal); CHKERRA(ierr);

      argl.clear();
      state.set_time(time);
      state_old.set_time(time_old);
      argl.arg_add(&state,IN_VECTOR|USE_TIME_DATA,"state");
      argl.arg_add(&state_old,IN_VECTOR|USE_TIME_DATA,"state_old");
      argl.arg_add(&res,OUT_VECTOR,"res");
      argl.arg_add(&res_delta,OUT_VECTOR,"res_delta");
      argl.arg_add(Ap,OUT_MATRIX|PFMAT,"A");
      argl.arg_add(&hmin,VECTOR_MIN,"hmin");
      argl.arg_add(&glob_param,USER_DATA,"glob_param");

      debug.trace("Before residual computation...");
      ierr = assemble(mesh,argl,dofmap,"comp_mat_res",&time);
      CHKERRA(ierr);
      debug.trace("After residual computation.");

      MPI_Allreduce(&mmv_functional, &mmv_functional_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&mmv_dfd, &mmv_dfd_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&mmv_d2fd, &mmv_d2fd_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&min_quality, &min_quality_g, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&min_volume, &min_volume_g, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      // tangled_mesh   = (min_volume > 0. ? 0 : 1);
      tangled_mesh   = (mmv_delta <= 1.e-10 ? 0 : 1);
      local_mmv_dfd  = mmv_dfd_g;
      local_mmv_d2fd = mmv_d2fd_g;

      if (!print_linear_system_and_stop || solve_system) {
	debug.trace("Before solving linear system...");
	ierr = Ap->solve(res,dx); CHKERRA(ierr); 
	debug.trace("After solving linear system.");
      }

      if (print_linear_system_and_stop) {
	PetscPrintf(PETSC_COMM_WORLD,
		    "Printing residual and matrix for"
		    " debugging and stopping.\n");
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,
				    "system.dat",&matlab); CHKERRA(ierr);
	ierr = PetscViewerSetFormat_WRAPPER(matlab,
					    PETSC_VIEWER_ASCII_MATLAB,
					    "atet"); CHKERRA(ierr);
	ierr = Ap->view(matlab);
	ierr = PetscViewerSetFormat_WRAPPER(matlab,
					    PETSC_VIEWER_ASCII_MATLAB,
					    "res"); CHKERRA(ierr);
	ierr = VecView(res,matlab);

	if (solve_system) {
	  ierr = PetscViewerSetFormat_WRAPPER(matlab,
					      PETSC_VIEWER_ASCII_MATLAB,"dx");
	  CHKERRA(ierr);
	  ierr = VecView(dx,matlab);
	}
	
	PetscFinalize();
	exit(0);

      }	

      double normres;
      ierr  = VecNorm(res,NORM_2,&normres); CHKERRA(ierr);
      if (check_for_inf) {
        PETSCFEM_ASSERT0(isfinite(normres),
                         "Detected Inf or NaN values in residual vector");  
      }
      if (inwt==0) normres_external = normres;
      PetscPrintf(PETSC_COMM_WORLD,
		  "Newton subiter %d, norm_res  = %10.3e\n",
		  inwt,normres);

      // substep update
      double relfac;
      int inwt_cum=0,nrf_indx=1;
      while (1) {
	if ((unsigned int)nrf_indx >= newton_relaxation_factor.size()) break;
	inwt_cum += int(newton_relaxation_factor[nrf_indx]);
	if (inwt_cum > inwt) break;
	nrf_indx += 2;
      }
      relfac = newton_relaxation_factor[nrf_indx-1];
      if (relfac!=1.) PetscPrintf(PETSC_COMM_WORLD,
				  "relaxation factor %f\n",relfac);
      scal= relfac/alpha;
      ierr = VecAXPY(x,scal,dx);

      ierr = VecDot(res_delta,dx,&resdelta_dot_dx);  CHKERRA(ierr);

      int    counter=0;
      if (use_line_search_mmv){
	double mmv_functional_ini = mmv_functional_g;
	double dif_functional,dx_dot_res;
	ierr = VecDot(res,dx,&dx_dot_res);  CHKERRA(ierr);
	tarea = 0;
	relfac = 1.;
	while (1){
	  mmv_functional = mmv_dfd = mmv_d2fd = 0.;
	  min_quality = 1.;
	  min_volume  = 1.e10;
	  argl.clear();
	  argl.arg_add(&state,IN_VECTOR|USE_TIME_DATA,"state");
	  argl.arg_add(&state_old,IN_VECTOR|USE_TIME_DATA,"state_old");
	  argl.arg_add(&resp,OUT_VECTOR,"res");
	  argl.arg_add(&res_delta,OUT_VECTOR,"res_delta");
	  argl.arg_add(Ap,OUT_MATRIX|PFMAT,"A");
	  argl.arg_add(&hmin,VECTOR_MIN,"hmin");
	  argl.arg_add(&glob_param,USER_DATA,"glob_param");
	
	  scal=0;
	  ierr = VecSet(resp,scal); CHKERRA(ierr);
	  ierr = Ap->clean_mat(); CHKERRA(ierr); 
	  debug.trace("Before residual computation...");
	  ierr = assemble(mesh,argl,dofmap,"comp_mat_res",&time);
	  CHKERRA(ierr);
	  debug.trace("After residual computation.");

	  MPI_Allreduce(&mmv_functional, &mmv_functional_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  MPI_Allreduce(&mmv_dfd, &mmv_dfd_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  MPI_Allreduce(&mmv_d2fd, &mmv_d2fd_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  MPI_Allreduce(&min_quality, &min_quality_g, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	  MPI_Allreduce(&min_volume, &min_volume_g, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	  if (min_volume_g > 0.) break;
	  dif_functional = fabs(mmv_functional_g)-
	    fabs(mmv_functional_ini+relfac*ls_eta_parameter*dx_dot_res);
	  if (dif_functional < 0.) break;

	  relfac /= ls_relfac_fraction;
	  scal = (1.-ls_relfac_fraction)/pow(ls_relfac_fraction,++counter);
	  ierr = VecAXPY(x,scal,dx);  CHKERRA(ierr);

	  PetscPrintf(PETSC_COMM_WORLD,"relaxation factor %f\n",relfac);
	  if (counter > 10) break;
	}
      }

      if (relfac==1.){
	// Incremento de delta
	mmv_delta_increment = fabs((local_mmv_dfd-0.0*resdelta_dot_dx)/local_mmv_d2fd);
	if (mmv_delta_increment > (1.-delta_fraction)*mmv_delta) 
	  mmv_delta_increment = (1.-delta_fraction)*mmv_delta;
	mmv_delta += -1.*mmv_delta_increment;
	if (min_volume_g <= 0. && mmv_delta < 1e-12)
	  mmv_delta = 1e-12;
      }
      if (counter > 10) mmv_delta *= 1.2;

#if 1
      PetscPrintf(PETSC_COMM_WORLD,"delta: %12.10f\n",mmv_delta);
      PetscPrintf(PETSC_COMM_WORLD,"Min quality: %f\n",min_quality_g);
      PetscPrintf(PETSC_COMM_WORLD,"Min quality old: %f\n",min_quality_old);
#endif
      
      if (fabs(min_quality_g-min_quality_old)/min_quality_g < tol_rel_quality 
      	  && min_quality_g > 0.) break;

      min_quality_old = min_quality_g;

#if 0
      ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD); CHKERRA(ierr);
      PetscFinalize();
      exit(0);
#endif

      if (normres < tol_newton) break;

      // fixme:= SHOULD WE CHECK HERE FOR NEWTON CONVERGENCE?

    } // end of loop over Newton subiteration (inwt)

    if (normres_external < tol_newton) {
      PetscPrintf(PETSC_COMM_WORLD,
		  "Tolerance on newton loop reached:  "
		  "|| R ||_0,  norm_res =%g < tol = %g\n",
		  normres_external,tol_newton);
      break;
    }

    // error difference
    scal = -1.0;
    ierr = VecAXPY(dx_step,scal,x);
    ierr  = VecNorm(dx_step,NORM_2,&norm); CHKERRA(ierr);
    if (check_for_inf) {
      PETSCFEM_ASSERT0(isfinite(norm),
                       "Detected Inf or NaN values in state vector");  
    }
    PetscPrintf(PETSC_COMM_WORLD,"============= delta_u = %10.3e\n",norm);
    print_vector_rota(save_file_pattern.c_str(),x,dofmap,&time,
		      tstep-1,nsaverot,nrec,nfile);
  
    if (nsave && tstep % nsave == 0) {
      print_vector(save_file.c_str(),x,dofmap,&time);
      if (print_residual) 
	print_vector(save_file_res.c_str(),res,dofmap,&time);
    }


    // Compute gathered quantities, for instance total force on walls
    if (ngather>0) {
      gather_values.resize(ngather,0.);
      for (int j=0; j<ngather; j++) gather_values[j] = 0.;
      arglf.clear();
      arglf.arg_add(&state,IN_VECTOR|USE_TIME_DATA);
      arglf.arg_add(&state_old,IN_VECTOR|USE_TIME_DATA);
      arglf.arg_add(&gather_values,VECTOR_ADD);
      ierr = assemble(mesh,arglf,dofmap,"gather",&time);
      CHKERRA(ierr);
    }

    hook_list.time_step_post(time.time(),tstep,gather_values);

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

    if (print_some_file!="" && tstep % nsome == 0) {
      print_some(save_file_some.c_str(),x,dofmap,node_list,&time);
    }

  }

  hook_list.close();

  if (report_option_access && MY_RANK==0) TextHashTable::print_stat();
  print_vector(save_file.c_str(),x,dofmap,&time);

  // ierr = VecDestroy(&x); CHKERRA(ierr); 
  ierr = VecDestroy(&x); CHKERRA(ierr); 
  ierr = VecDestroy(&xold); CHKERRA(ierr); 
  ierr = VecDestroy(&dx); CHKERRA(ierr); 
  ierr = VecDestroy(&dx_step); CHKERRA(ierr); 
  ierr = VecDestroy(&res); CHKERRA(ierr); 
  ierr = VecDestroy(&res_delta); CHKERRA(ierr); 

  DELETE_SCLR(Ap);
  DELETE_SCLR(dofmap);
  DELETE_SCLR(mesh);
#ifdef DEBUG_MALLOC_USE
  fclose(malloc_log);
#endif
  PetscFinalize();
  exit(0);
}
