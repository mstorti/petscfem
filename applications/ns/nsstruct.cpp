//__INSERT_LICENSE__
//$Id new structure module with displacement formulation - rodrigop -0300$
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
#define PetscViewerSetFormat_WRAPPER(viewer,format,name)	\
  PetscViewerSetFormat(viewer,format)

#include <applications/ns/nsi_tet.h>
static char help[] = "PETSc-FEM Navier Stokes module\n\n";

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Creates hooks depending on the name. 
    @param name (input) the name of the hook. 
    @return a pointer to the hook. */ 
Hook *ns_hook_factory(const char *name);

extern vector<double> data_pts;
extern vector<ElemToPtr> elemset_pointer;

//debug:=
extern int TSTEP;
extern int fractional_step;
extern int reuse_mat;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "struct_main"
int struct_main() {

  Vec x, xph, xmh, xold, dx, dx_step, res; // approx solution, RHS, residual
  PetscViewer matlab;
  PFMat *A_tet=NULL, *A_tet_c=NULL, *A_mom=NULL, 
    *A_poi=NULL, *A_prj=NULL;	// linear system matrix 
  double  norm, *sol, scal;	// norm of solution error
  int     ierr, i, n = 10, size, node,
    jdof, k, kk, nfixa,
    kdof, ldof, lloc, nel, nen, neq, nu,
    myrank;
  PetscTruth flg;
  // Initialize time
  // time = t^{n+1} -> to be computed
  // time_old = t^{n} (known)
  // time_ph = t^{n+1/2} (known)
  // time_mh = t^{n-1/2} (known)
  // In fact (xold, xmh, xph) are a way yo rewrite
  // (the values of x, \dot x, and \ddot x) at t^{n}, i.e.
  // \dot x = (xph-xmh)/Dt
  // \ddot x = (xph-2*xold+xmh)/Dt^2
  Time time,time_old,time_ph,time_mh,time_star; 
  GlobParam glob_param;
  GLOB_PARAM = &glob_param;
  string save_file_res;
  BasicObject_application_factory = &BasicObject_ns_factory;
  
  // ierr = MatCreateShell(PETSCFEM_COMM_WORLD,int m,int n,int M,int N,void *ctx,Mat *A)
  char fcase[FLEN+1],output_file[FLEN+1];
  Dofmap *dofmap;
  Mesh *mesh;
  // arglf:= argument list for computing gathered quantities as forces
  arg_list argl, arglf;
  vector<double> hmin;
  hmin.resize(1);

  print_copyright();
  PetscPrintf(PETSCFEM_COMM_WORLD,"-------- Structure Code for Fluid "
              "Structure Interaction module ---------\n");

  Debug debug(0,PETSCFEM_COMM_WORLD);
  GLOBAL_DEBUG = &debug;

  //  if (size != 1) SETERRA(1,0,"This is a uniprocessor example only!");
  ierr = PetscOptionsGetString(PETSC_NULL,"-case",fcase,FLEN,&flg);
  CHKERRA(ierr);
  if (!flg) {
    PetscPrintf(PETSCFEM_COMM_WORLD,
		"Option \"-case <filename>\" not passed to PETSc-FEM!!\n");
    PetscFinalize();
    exit(0);
  }

  ierr = PetscOptionsGetString(PETSC_NULL,"-o",output_file,FLEN,&flg);
  CHKERRA(ierr);
  if (flg) { 
    PetscPrintf(PETSCFEM_COMM_WORLD,"PETSc-FEM: NS module: "
		"redirecting output to \"%s\"\n",output_file);
    FILE *new_stdout = fopen(output_file,"w");
    if (!new_stdout) {
      PetscPrintf(PETSCFEM_COMM_WORLD,"error redirecting output. "
		  "Couldn't open \"%s\"\n",output_file);
    } else {
      fclose(stdout);
      stdout = new_stdout;
    }
  }

  // Read the mesh
  ierr = read_mesh(mesh,fcase,dofmap,neq,SIZE,MY_RANK); CHKERRA(ierr); 
  PetscPrintf(PETSCFEM_COMM_WORLD,"After readmesh...\n");

  GLOBAL_OPTIONS = mesh->global_options;

  //o Time to start computations
  TGETOPTDEF(GLOBAL_OPTIONS,double,start_comp_time,0.);
  time.set(start_comp_time);
  time_mh.set(start_comp_time);
  time_ph.set(start_comp_time);
  State 
    state(x,time), 
    state_old(xold,time_old),
    state_mh(xmh,time_mh),
    state_ph(xph,time_ph);

#if 0
  //o If set, redirect output to this file.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,stdout_file,none);
  if (!strcmp(stdout_file.c_str(),"")) {
    fclose(stdout);
    stdout = fopen(stdout_file.c_str(),"w");
  }
#endif

#if 0
  // Check that `use_displacement_formulation' was not set by user
  TGETOPTDEF(GLOBAL_OPTIONS,int,use_displacement_formulation,-1);
  PETSCFEM_ASSERT0(use_displacement_formulation==-1,
                   "use_displacement_formulation must not be set by user");  
  // Use new formulation 
  GLOBAL_OPTIONS->add_entry("use_displacement_formulation","1",0);
#endif
  GLOBAL_OPTIONS->add_entry("use_arg_handles","1",0);

  //o The number of outer stages for convergence in coupled problems.
  TGETOPTDEF(GLOBAL_OPTIONS,int,nstage,1);
  int fractional_step = 0;
  //o Reuse matrices
  TGETOPTDEF_ND(GLOBAL_OPTIONS,int,reuse_mat,0);
  //o Fractional step uses symmetric matrices (only CG iterative KSP).
  TGETOPTDEF(GLOBAL_OPTIONS,int,fractional_step_use_petsc_symm,1);
  //o Solver combination for the fractional step method. May be #iisd#, 
  //  #lu#, #global_gmres#. 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,fractional_step_solver_combo,"iisd");

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

  HookList hook_list;
  hook_list.init(*mesh,*dofmap,ns_hook_factory);

  //o Dimension of the problem.
  GETOPTDEF(int,ndim,3);

  //o If 0: compute #wall_data# info only once. Otherwise
  //   refresh each #update_wall_data# steps. 
  GETOPTDEF(int,update_wall_data,0);

  //o Scales displacement for ALE-like mesh relocation. 
  GETOPTDEF(double,displ_factor,0.1);
  //o Number of inner iterations for the global non-linear
  // Newton  problem. 
  GETOPTDEF(int,nnwt,1);
  //o Update jacobian each $n$-th time step. 
  GETOPTDEF(int,update_jacobian_steps,0);
#define INF INT_MAX
  //o Update jacobian each $n$-th time step. 
  GETOPTDEF(int,update_jacobian_start_steps,INF);
#undef INF
  //o Tolerance to solve the non-linear system (global Newton).
  GETOPTDEF(double,tol_newton,1e-8);

#define INF INT_MAX
  //o Update jacobian only until n-th Newton subiteration. 
  // Don't update if null. 
  GETOPTDEF(int,update_jacobian_iters,1);
  assert(update_jacobian_iters>=1);
  //o Update jacobian each $n$-th Newton iteration
  GETOPTDEF(int,update_jacobian_start_iters,INF);
  assert(update_jacobian_start_iters>=0);
#undef INF

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

  //o Use the LES/Smagorinsky turbulence model. 
  GETOPTDEF(int,LES,0);
  //o If {\tt A\_van\_Driest=0} then the van Driest
  //    damping factor is not used 
  GETOPTDEF(int,A_van_Driest,0);

  if(A_van_Driest>0) { 
    PetscPrintf(PETSCFEM_COMM_WORLD,"--- Don't forget to refresh Wall_Data -- \n");
    PetscPrintf(PETSCFEM_COMM_WORLD,"--- using update_wall_data global option -- \n");
  }

  //o frequency of Poisson matrix recomputation
  GETOPTDEF(int,freq_update_mat_poi,1000000);

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
  //o The pattern to generate the file name to save in for Newton iterations.
  //  If not set, Newton iterations are not saved.
  // the rotary save mechanism.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file_pattern_nwt,<none>);
  int save_newton_iters = (save_file_pattern_nwt!="<none>");

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

  LPFilterGroup filter(GLOBAL_OPTIONS,sx,Dt);
  
  // Use IISD (Interface Iterative Subdomain Direct) or not.
  // A_tet = (use_iisd ? &IISD_A_tet : &PETSc_A_tet);
  A_tet = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());
  A_tet_c = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());

  ierr = VecDuplicate(x,&xold); CHKERRA(ierr);
  ierr = VecDuplicate(x,&xmh); CHKERRA(ierr);
  ierr = VecDuplicate(x,&xph); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx_step); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx); CHKERRA(ierr);
  ierr = VecDuplicate(x,&res); CHKERRA(ierr);

  // Just for debugging
  ierr = VecSet(x,Dt); CHKERRA(ierr);
  ierr = VecSet(xold,0.0); CHKERRA(ierr);
  ierr = VecSet(xmh,-Dt/2.0); CHKERRA(ierr);
  ierr = VecSet(xph,+Dt/2.0); CHKERRA(ierr);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // initialize state vectors
  scal=0;
  ierr = VecSet(x,scal); CHKERRA(ierr);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Compute  profiles
  debug.trace("Computing profiles...");
  argl.clear();
  argl.arg_add(A_tet,PROFILE|PFMAT,"A");
  ierr = assemble(mesh,argl,dofmap,"comp_prof",&time); CHKERRA(ierr); 
  debug.trace("After computing profile.");
  int update_jacobian_step=0;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  debug.trace("Before reading initial vector...");
  ierr = opt_read_vector(mesh,x,dofmap,MY_RANK); CHKERRA(ierr);
  debug.trace("After reading initial vector...");

  // update_jacobian_this_step:= Flags whether this step the
  // jacobian should be updated or not 
  int update_jacobian_this_step,
    update_jacobian_this_iter, tstep_start=1;
  for (int tstep=tstep_start; tstep<=nstep; tstep++) {
    TSTEP=tstep; //debug:=
    time_old.set(time.time());
    time_mh.set(time.time());
    time_mh.inc(-Dt/2);
    time_ph.set(time.time());
    time_ph.inc(+Dt/2);
    time_star.set(time.time()+alpha*Dt);
    time.inc(Dt);
    if (!MY_RANK) printf("Time step: %d, time: %g %s\n",
			 tstep,time.time(),(steady ? " (steady) " : ""));

    hook_list.time_step_pre(time_star.time(),tstep);
    ierr = VecCopy(x,xold);
    
    for (int stage=0; stage<nstage; stage++) {
      
      PetscPrintf(PETSCFEM_COMM_WORLD,
		  " --------------------------------------\n"
		  "Stage in Navier-Stokes #: %d / %d \n"
		  " --------------------------------------\n",
		  stage,nstage);
      
      hook_list.stage("stage_pre",stage,time_star.time(),&gather_values);

      // Jacobian update logic
      update_jacobian_this_step = (tstep < update_jacobian_start_steps) 
	|| ((tstep-update_jacobian_start_steps) % update_jacobian_steps == 0);
    
      // Inicializacion del paso
      // if (stage>0) ierr = VecCopy(xold,x);
      ierr = VecCopy(x,dx_step);

      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
      // NON LINEAR ALGORITHM
      double normres_external=NAN;
      for (int inwt=0; inwt<nnwt; inwt++) {

        glob_param.inwt = inwt;
        // Initialize step

        // update_jacobian_this_iter:= flags whether the Jacobian is
        // factored in this iter or not Jacobian update logic
        update_jacobian_this_iter = update_jacobian_this_step &&
          ( (inwt < update_jacobian_start_iters) 
            || ((inwt - update_jacobian_start_iters) % update_jacobian_iters == 0) );

        if (update_jacobian_this_iter) {
          ierr = A_tet->clean_mat(); CHKERRA(ierr); 
        }

        scal=0;
        ierr = VecSet(res,scal); CHKERRA(ierr);
        if (update_jacobian_this_iter) {
          // ierr = A_tet->clean_mat(); CHKERRA(ierr); 
        }

        argl.clear();
        state.set_time(time);
        state_old.set_time(time_old);
        argl.arg_add(&state,IN_VECTOR|USE_TIME_DATA,"state");
        argl.arg_add(&state_old,IN_VECTOR|USE_TIME_DATA,"state_old");
        argl.arg_add(&state_mh,IN_VECTOR|USE_TIME_DATA,"state_mh");
        argl.arg_add(&state_ph,IN_VECTOR|USE_TIME_DATA,"state_ph");
        argl.arg_add(&res,OUT_VECTOR,"res");
        if (update_jacobian_this_iter) argl.arg_add(A_tet,OUT_MATRIX|PFMAT,"A");
        argl.arg_add(&hmin,VECTOR_MIN,"hmin");
        argl.arg_add(&glob_param,USER_DATA,"glob_param");

        const char *jobinfo = (update_jacobian_this_iter ? "comp_mat_res" : "comp_res");

        // In order to measure performance
        if (measure_performance) {
          ierr = measure_performance_fun(mesh,argl,dofmap,jobinfo,
                                         &time_star); CHKERRA(ierr);
          PetscFinalize();
          exit(0);
        }

        debug.trace("Before residual computation...");
        ierr = assemble(mesh,argl,dofmap,jobinfo,&time_star);
        CHKERRA(ierr);
        debug.trace("After residual computation.");

#if 0
        ierr = PetscViewerASCIIOpen(PETSCFEM_COMM_WORLD,
                                    "system.dat",&matlab); CHKERRA(ierr);
        ierr = PetscViewerSetFormat(matlab,
                                    PETSC_VIEWER_ASCII_MATLAB); CHKERRA(ierr);
        // ierr = PetscViewerSetFilename(matlab,"a_tet"); CHKERRA(ierr);
        ierr = A_tet->view(matlab); CHKERRA(ierr); 
#endif

#if 0 //debug:=
        ierr = PetscViewerASCIIOpen(PETSCFEM_COMM_WORLD,
                                    "system.dat",&matlab); CHKERRA(ierr);
        ierr = PetscViewerSetFormat_WRAPPER(matlab,
                                            PETSC_VIEWER_ASCII_MATLAB,"x"); CHKERRA(ierr);
        ierr = VecView(x,matlab);
        ierr = PetscViewerSetFormat_WRAPPER(matlab,
                                            PETSC_VIEWER_ASCII_MATLAB,"res"); CHKERRA(ierr);
        ierr = VecView(res,matlab);
        PetscFinalize();
        exit(0);
#endif

        PetscViewer matlab;
        if (verify_jacobian_with_numerical_one) {
          PETSCFEM_ERROR0("NOT IMPLEMENTED YET");  
          PETSCFEM_ERROR0("Not implemented for this version of struct code");  
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
          argl.arg_add(&x,PERT_VECTOR,"state");
          argl.arg_add(&xold,IN_VECTOR,"state_old");
          argl.arg_add(A_tet_c,OUT_MATRIX_FDJ|PFMAT,"A");

          update_jacobian_step++;
          if (update_jacobian_step >= update_jacobian_steps) 
            update_jacobian_step =0;
          if (update_jacobian_this_iter) argl.arg_add(A_tet,OUT_MATRIX|PFMAT);
          argl.arg_add(&hmin,VECTOR_MIN);
          argl.arg_add(&glob_param,USER_DATA);
          ierr = assemble(mesh,argl,dofmap,jobinfo,
                          &time_star); CHKERRA(ierr);

          ierr = PetscViewerSetFormat_WRAPPER(matlab,
                                              PETSC_VIEWER_ASCII_MATLAB,"atet_fdj"); CHKERRA(ierr);
          ierr = A_tet_c->view(matlab); CHKERRQ(ierr); 

          PetscFinalize();
          exit(0);
        }

        if (!print_linear_system_and_stop || solve_system) {
          debug.trace("Before solving linear system...");
          ierr = A_tet->solve(res,dx); CHKERRA(ierr); 
          debug.trace("After solving linear system.");
        }

        if (print_linear_system_and_stop) {
          PetscPrintf(PETSCFEM_COMM_WORLD,
                      "Printing residual and matrix for"
                      " debugging and stopping.\n");
          ierr = PetscViewerASCIIOpen(PETSCFEM_COMM_WORLD,
                                      "system.dat",&matlab); CHKERRA(ierr);
          ierr = PetscViewerSetFormat_WRAPPER(matlab,
                                              PETSC_VIEWER_ASCII_MATLAB,"atet"); CHKERRA(ierr);
          ierr = A_tet->view(matlab);
          ierr = PetscViewerSetFormat_WRAPPER(matlab,
                                              PETSC_VIEWER_ASCII_MATLAB,"res"); CHKERRA(ierr);
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
        if (inwt==0) normres_external = normres;
        PetscPrintf(PETSCFEM_COMM_WORLD,
                    "Newton subiter %d, norm_res  = %10.3e, update Jac. %d\n",
                    inwt,normres,update_jacobian_this_iter);

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
        if (relfac!=1.) PetscPrintf(PETSCFEM_COMM_WORLD,
                                    "relaxation factor %f\n",relfac);
        scal= relfac/alpha;
        ierr = VecAXPY(x,scal,dx); CHKERRA(ierr);

#if 0
        ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD); CHKERRA(ierr);
        PetscFinalize();
        exit(0);
#endif

        if (save_newton_iters) {
          char save_file_nwt[200];
          sprintf(save_file_nwt,save_file_pattern_nwt.c_str(),
                  tstep,inwt);
	    
          print_vector(save_file_nwt,x,dofmap,&time);
        }

        // fixme:= SHOULD WE CHECK HERE FOR NEWTON CONVERGENCE?	
        if (normres_external < tol_newton) {
          PetscPrintf(PETSCFEM_COMM_WORLD,
                      "Tolerance on newton loop reached:  "
                      "|| R ||_0,  norm_res =%g < tol = %g\n",
                      normres_external,tol_newton);
          break;
        }	
      } // end of loop over Newton subiteration (inwt)
    
      // error difference
      scal = -1.0;
      ierr = VecAXPY(dx_step,scal,x);
      ierr  = VecNorm(dx_step,NORM_2,&norm); CHKERRA(ierr);
      PetscPrintf(PETSCFEM_COMM_WORLD,"============= delta_u = %10.3e\n",norm);
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
	ierr = assemble(mesh,arglf,dofmap,"gather",&time_star);
	CHKERRA(ierr);
      }

      int converged = 0;
      hook_list.stage("stage_post",stage,
		      time_star.time(),&converged);
      if (converged) break;
      
    } // end for stage
    

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

    filter.update(time);
    if (print_some_file!="" && tstep % nsome == 0) {
      print_some(save_file_some.c_str(),x,dofmap,node_list,&time);
      filter.print_some("filter.some.tmp",dofmap,node_list);
    }

  }

  hook_list.close();

  if (report_option_access && MY_RANK==0) TextHashTable::print_stat();
  print_vector(save_file.c_str(),x,dofmap,&time);

  // ierr = VecDestroy(x); CHKERRA(ierr); 
  ierr = VecDestroy(x); CHKERRA(ierr); 
  ierr = VecDestroy(xold); CHKERRA(ierr); 
  ierr = VecDestroy(dx); CHKERRA(ierr); 
  ierr = VecDestroy(dx_step); CHKERRA(ierr); 
  ierr = VecDestroy(res); CHKERRA(ierr); 

  delete A_tet;
  delete A_tet_c;

  DELETE_SCLR(dofmap);
  DELETE_SCLR(mesh);
#ifdef DEBUG_MALLOC_USE
  fclose(malloc_log);
#endif
  PetscFinalize();
  exit(0);
}
