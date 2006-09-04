//__INSERT_LICENSE__
//$Id: ns.cpp,v 1.188 2006/09/04 00:21:00 mstorti Exp $
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
#include <src/iisdmatstat.h>

// PETSc now doesn't have the string argument that represents the variable name
// so that I will use this wrapper until I find how to set names in Ascii matlab viewers.
#define PetscViewerSetFormat_WRAPPER(viewer,format,name) \
          PetscViewerSetFormat(viewer,format)

#include <applications/ns/nsi_tet.h>
static char help[] = "PETSc-FEM Navier Stokes module\n\n";

extern int MY_RANK,SIZE;
WallData wall_data;

int fsi_main();
int struct_main();
int mmove_main();

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Creates hooks depending on the name. 
    @param name (input) the name of the hook. 
    @return a pointer to the hook. */ 
Hook *ns_hook_factory(const char *name);

vector<double> data_pts;
vector<ElemToPtr> elemset_pointer;

//debug:=
int TSTEP=0;
int fractional_step;
int reuse_mat;

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
#undef __FUNC__
GlobParam *GLOB_PARAM;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define __FUNC__ "main"
int main(int argc,char **args) {

  Vec x, xp, dx, xold, dx_step, res; // approx solution, RHS, residual
  PetscViewer matlab;
  PFMat *A_tet, *A_tet_c, *A_mom, *A_poi, *A_prj;;	// linear system matrix 
  double  norm, *sol, scal;	// norm of solution error
  int     ierr, i, n = 10, size, node,
    jdof, k, kk, nfixa,
    kdof, ldof, lloc, nel, nen, neq, nu,
    myrank;
  PetscTruth flg;
  // Initialize time
  Time time,time_old,time_star; 
  GlobParam glob_param;
  GLOB_PARAM = &glob_param;
  string save_file_res;
  BasicObject_application_factory = &BasicObject_ns_factory;
  
  // ierr = MatCreateShell(PETSC_COMM_WORLD,int m,int n,int M,int N,void *ctx,Mat *A)
  char fcase[FLEN+1],output_file[FLEN+1];
  Dofmap *dofmap;
  Mesh *mesh;
  // arglf:= argument list for computing gathered quantities as forces
  arg_list argl, arglf;
  vector<double> hmin;
  hmin.resize(1);

  PetscInitialize(&argc,&args,(char *)0,help);

  // Get MPI info
  MPI_Comm_size(PETSC_COMM_WORLD,&SIZE);
  MPI_Comm_rank(PETSC_COMM_WORLD,&MY_RANK);

#define CNLEN 100
  char code_name[CNLEN];
  ierr = PetscOptionsGetString(PETSC_NULL,"-code",code_name,CNLEN,&flg);
  printf("flg %d, code %s\n",flg,code_name);

  if (flg) {
    if (!strcmp(code_name,"fsi")) return fsi_main();
    if (!strcmp(code_name,"struct")) return struct_main();
    if (!strcmp(code_name,"mmove")) return mmove_main();
  }

  print_copyright();
  PetscPrintf(PETSC_COMM_WORLD,"-------- Navier-Stokes module ---------\n");

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
  State state(x,time),statep(xp,time),state_old(xold,time_old);

#if 0
  //o If set, redirect output to this file.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,stdout_file,none);
  if (!strcmp(stdout_file.c_str(),"")) {
    fclose(stdout);
    stdout = fopen(stdout_file.c_str(),"w");
  }
#endif

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
  for (int j=1; j<newton_relaxation_factor.size(); j += 2) {
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

  //o Sets the save frequency in iterations 
  GETOPTDEF(int,nsave,10);
  //o Sets the frequency save for the ``rotary save'' mechanism. 
  //i_tex ../../doc/nsdoc.tex rotary_save
  GETOPTDEF(int,nsaverot,100);
  //o Sets the number of states saved in a given file
  // in the ``rotary save'' mechanism (see \ref{sec:rotary_save}
  GETOPTDEF(int,nrec,1000000);
  //o Sets the number of files in the ``rotary save'' mechanism. 
  // (see \ref{sec:rotary_save})
  GETOPTDEF(int,nfile,1);

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
    fprintf(gather_file_f,"");
    fclose(gather_file_f);
  }

  //o Use the LES/Smagorinsky turbulence model. 
  GETOPTDEF(int,LES,0);
  //o If {\tt A\_van\_Driest=0} then the van Driest
  //    damping factor is not used 
  GETOPTDEF(int,A_van_Driest,0);

  if(A_van_Driest>0) { 
    PetscPrintf(PETSC_COMM_WORLD,"--- Don forget to refresh Wall_Data -- \n");
    PetscPrintf(PETSC_COMM_WORLD,"--- using update_wall_data global option -- \n");
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

  LPFilterGroup filter(GLOBAL_OPTIONS,sx,Dt);
  
  // Use IISD (Interface Iterative Subdomain Direct) or not.
  // A_tet = (use_iisd ? &IISD_A_tet : &PETSc_A_tet);
  if (!fractional_step) {
    A_tet = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());
    A_tet_c = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());
  } else {
    if (!fractional_step_use_petsc_symm) {
      A_mom = PFMat::dispatch(dofmap->neq,*dofmap,solver_mom.c_str());
      A_poi = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());
      A_poi->set_option("KSP_method","cg");
      A_poi->set_option("preco_side","left");
      A_prj = PFMat::dispatch(dofmap->neq,*dofmap,solver_mom.c_str());
    } else {
      if (fractional_step_solver_combo=="iisd") {
	// IISD (Domain decomposition) iteration
	A_mom = PFMat::dispatch(dofmap->neq,*dofmap,"iisd");
	A_mom->set_option("preco_type","jacobi");
	A_mom->set_option("print_internal_loop_conv","1");
#define IISDMAT_PRINT_STATISTICS "0"
	A_mom->set_option("iisdmat_print_statistics",IISDMAT_PRINT_STATISTICS);
	A_mom->set_option("use_interface_full_preco_nlay","1");
	A_mom->set_option("block_uploading","2");

	A_poi = PFMat::dispatch(dofmap->neq,*dofmap,"iisd");
	A_poi->set_option("preco_type","jacobi");
	A_poi->set_option("print_internal_loop_conv","1");
	A_poi->set_option("block_uploading","0");
	A_poi->set_option("iisdmat_print_statistics",IISDMAT_PRINT_STATISTICS);
	A_poi->set_option("use_interface_full_preco_nlay","1");

	A_prj = PFMat::dispatch(dofmap->neq,*dofmap,"petsc");
	A_prj->set_option("preco_type","jacobi");
	A_prj->set_option("print_internal_loop_conv","1");
	A_prj->set_option("block_uploading","0");

      } else if (fractional_step_solver_combo=="global_gmres") {

	// Global PETSc iteration
	// A_mom->set_option("KSP_method",KSPBCGS);
	A_mom = PFMat::dispatch(dofmap->neq,*dofmap,"petsc");
	A_mom->set_option("KSP_method",KSPGMRES);
	A_mom->set_option("preco_side","left");
	A_mom->set_option("preco_type","jacobi");

#if 1
	A_poi = PFMat::dispatch(dofmap->neq,*dofmap,"petsc");
	A_poi->set_option("KSP_method",KSPCG);
	A_poi->set_option("preco_type","jacobi");
#else
	A_poi = PFMat::dispatch(dofmap->neq,*dofmap,"iisd");
	A_poi->set_option("preco_type","jacobi");
	A_poi->set_option("print_internal_loop_conv","1");
	A_poi->set_option("block_uploading","0");
	A_poi->set_option("iisdmat_print_statistics",IISDMAT_PRINT_STATISTICS);
	A_poi->set_option("use_interface_full_preco_nlay","1");
#endif

	A_prj = PFMat::dispatch(dofmap->neq,*dofmap,"petsc_symm");
	A_prj->set_option("KSP_method",KSPCG);
	A_prj->set_option("preco_type","jacobi");
	A_prj->set_option("block_uploading","0");

      } else if (fractional_step_solver_combo=="lu") {
#define SET_SOLVER_OPTIONS(solv)				\
    solv = PFMat::dispatch(dofmap->neq,*dofmap,"petsc");	\
      solv->set_option("preco_type","lu");			\
      solv->set_option("block_uploading","0");			\
      solv->set_option("KSP_method",KSPRICHARDSON); 

      SET_SOLVER_OPTIONS(A_mom);
      SET_SOLVER_OPTIONS(A_poi);
      SET_SOLVER_OPTIONS(A_prj);
      } else {
	PETSCFEM_ERROR("Unknown fractional step solver %d\n",
		       fractional_step_solver_combo.c_str());
      }
    }
    ierr = VecDuplicate(x,&xp); CHKERRA(ierr);
  }

  ierr = VecDuplicate(x,&xold); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx_step); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx); CHKERRA(ierr);
  ierr = VecDuplicate(x,&res); CHKERRA(ierr);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // initialize state vectors
  scal=0;
  ierr = VecSet(&scal,x); CHKERRA(ierr);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Compute  profiles
  debug.trace("Computing profiles...");
  if (!fractional_step) {
    argl.clear();
    argl.arg_add(A_tet,PROFILE|PFMAT);
    ierr = assemble(mesh,argl,dofmap,"comp_mat",&time); CHKERRA(ierr); 

  } else {
    argl.clear();
    argl.arg_add(A_mom,PROFILE|PFMAT);
    argl.arg_add(A_poi,PROFILE|PFMAT);
    argl.arg_add(A_prj,PROFILE|PFMAT);
    ierr = assemble(mesh,argl,dofmap,"comp_mat_prof",&time); CHKERRA(ierr); 
  }
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
    time_star.set(time.time()+alpha*Dt);
    time.inc(Dt);
    if (!MY_RANK) printf("Time step: %d, time: %g %s\n",
			 tstep,time.time(),(steady ? " (steady) " : ""));


    //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    // Build octree for nearest neighbor calculation
    int refresh_wall_data = LES && A_van_Driest>0.;
    int rest_ok = update_wall_data 
      && ((tstep-tstep_start) % update_wall_data)==0;
    refresh_wall_data &= (tstep==tstep_start || rest_ok);
    if (refresh_wall_data) {
      PetscPrintf(PETSC_COMM_WORLD,"--- Refreshing Wall_Data -- \n");
#ifdef USE_ANN
      argl.clear();
      argl.arg_add(&data_pts,USER_DATA);
      argl.arg_add(&elemset_pointer,USER_DATA);
      ierr = assemble(mesh,argl,dofmap,"build_nneighbor_tree",&time); CHKERRQ(ierr); 
      PetscSynchronizedFlush(PETSC_COMM_WORLD);
      dvector<double> buff;
      buff.mono(data_pts.size());
      ierr = MPI_Allreduce(&*data_pts.begin(),buff.buff(),
			   data_pts.size(),MPI_DOUBLE,
			   MPI_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);
      for (int j=0; j<data_pts.size(); j++) data_pts[j] = buff.ref(j);
      buff.clear();

      PetscPrintf(PETSC_COMM_WORLD,"After nearest neighbor tree.\n");

      wall_data.init(&data_pts,&elemset_pointer,ndim);

      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
      // Find nearest neighbor for each volume element
      argl.clear();
      argl.arg_add(&wall_data,USER_DATA);
      ierr = assemble(mesh,argl,dofmap,"get_nearest_wall_element",
		      &time); CHKERRA(ierr); 
#else
      PETSCFEM_ERROR0("Not compiled with ANN library!!\n");
#endif
    }

#if 0 && defined USE_ANN
    ANNpoint point = annAllocPt(ndim);
    ANNidx nn_idx;
    ANNpoint nn = annAllocPt(ndim);
    ANNdist dist;
    int elem;


    while (1) {
      printf("\n\nEnter point: > ");
      scanf("%lf %lf %lf",&(point[0]),&(point[1]),&(point[2]));
      printf("\n");

      wall_data->nearest(point,elemset,elem,nn_idx,nn,dist);
      printf("Nearest neighbor id point %d, coords: %f %f %f, elem %d, elemset %p\n",
	     nn_idx,nn[0],nn[1],nn[2],elem,elemset);
    }
#endif

    hook_list.time_step_pre(time_star.time(),tstep);
    // Jacobian update logic
    update_jacobian_this_step = (tstep < update_jacobian_start_steps) 
      || ((tstep-update_jacobian_start_steps) % update_jacobian_steps == 0);
    
    // Inicializacion del paso
    ierr = VecCopy(x,dx_step);
    ierr = VecCopy(x,xold);
    
    if (!fractional_step) {
      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
      // TET ALGORITHM
      
      double normres_external;
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

	// Compute wall stresses
	argl.clear();
	argl.arg_add(&x,IN_VECTOR);
	argl.arg_add(&xold,IN_VECTOR);
	ierr = assemble(mesh,argl,dofmap,"comp_shear_vel",
			&time_star); CHKERRA(ierr);

	// Communicate wall stresses amoung different processors
	// This is obsolete. Now we control this from the loop
	// over elements in the `wall' elemset. 

	// Mon Oct 23 19:13:39 ART 2000. Now I think that I need this
	// because the at each processor there is not a global version
	// of the state.
	if (LES && A_van_Driest>0. && SIZE>1) {
	  argl.clear();
	  ierr = assemble(mesh,argl,dofmap,"communicate_shear_vel",
			  &time_star); CHKERRA(ierr);
	}

	scal=0;
	ierr = VecSet(&scal,res); CHKERRA(ierr);
	if (update_jacobian_this_iter) {
	  // ierr = A_tet->clean_mat(); CHKERRA(ierr); 
	}

	argl.clear();
	state.set_time(time);
	state_old.set_time(time_old);
	argl.arg_add(&state,IN_VECTOR|USE_TIME_DATA);
	argl.arg_add(&state_old,IN_VECTOR|USE_TIME_DATA);
	argl.arg_add(&res,OUT_VECTOR);
	if (update_jacobian_this_iter) argl.arg_add(A_tet,OUT_MATRIX|PFMAT);
	argl.arg_add(&hmin,VECTOR_MIN);
	argl.arg_add(&glob_param,USER_DATA);
	argl.arg_add(&wall_data,USER_DATA);

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
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,
			       "system.dat",&matlab); CHKERRA(ierr);
	ierr = PetscViewerSetFormat(matlab,
			       PETSC_VIEWER_ASCII_MATLAB); CHKERRA(ierr);
	// ierr = PetscViewerSetFilename(matlab,"a_tet"); CHKERRA(ierr);
	ierr = A_tet->view(matlab); CHKERRA(ierr); 
#endif

#if 0 //debug:=
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,
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
	  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,
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

	  update_jacobian_step++;
	  if (update_jacobian_step >= update_jacobian_steps) 
	    update_jacobian_step =0;
	  if (update_jacobian_this_iter) argl.arg_add(A_tet,OUT_MATRIX|PFMAT);
	  argl.arg_add(&hmin,VECTOR_MIN);
	  argl.arg_add(&glob_param,USER_DATA);
	  argl.arg_add(&wall_data,USER_DATA);
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
	  PetscPrintf(PETSC_COMM_WORLD,
		      "Printing residual and matrix for"
		      " debugging and stopping.\n");
	  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,
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
	PetscPrintf(PETSC_COMM_WORLD,
		    "Newton subiter %d, norm_res  = %10.3e, update Jac. %d\n",
		    inwt,normres,update_jacobian_this_iter);

	// update del subpaso
	double relfac;
	int inwt_cum=0,nrf_indx=1;
	while (1) {
	  if (nrf_indx >= newton_relaxation_factor.size()) break;
	  inwt_cum += int(newton_relaxation_factor[nrf_indx]);
	  if (inwt_cum > inwt) break;
	  nrf_indx += 2;
	}
	relfac = newton_relaxation_factor[nrf_indx-1];
	if (relfac!=1.) PetscPrintf(PETSC_COMM_WORLD,
				    "relaxation factor %f\n",relfac);
	scal= relfac/alpha;
	ierr = VecAXPY(&scal,dx,x);

#if 0
	ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD); CHKERRA(ierr);
	PetscFinalize();
	exit(0);
#endif

	// fixme:= SHOULD WE CHECK HERE FOR NEWTON CONVERGENCE?
	if (normres_external < tol_newton) {
	  PetscPrintf(PETSC_COMM_WORLD,
		      "Tolerance on newton loop reached:  "
		      "|| R ||_0,  norm_res =%g < tol = %g\n",
		      normres_external,tol_newton);
	  break;
	}

      } // end of loop over Newton subiteration (inwt)

    } else {
    
      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
      // FRACTIONAL STEP ALGORITHM
      // Inicializacion del paso
      ierr = VecCopy(x,dx_step);
      ierr = VecCopy(x,xold);

      // Compute wall stresses
      argl.clear();
      argl.arg_add(&x,IN_VECTOR);
      argl.arg_add(&xold,IN_VECTOR);
      ierr = assemble(mesh,argl,dofmap,"comp_shear_vel",
		      &time_star); CHKERRA(ierr);

      // Communicate wall stresses amoung different processors
      // This is obsolete. Now we control this from the loop
      // over elements in the `wall' elemset. 

	// Mon Oct 23 19:13:39 ART 2000. Now I think that I need this
	// because the at each processor there is not a global version
	// of the state.
      if (LES && A_van_Driest>0. && SIZE>1) {
	argl.clear();
	ierr = assemble(mesh,argl,dofmap,"communicate_shear_vel",
			&time_star); CHKERRA(ierr);
      }

#define PRE_STEP
#define POI_STEP
#define PRJ_STEP

#ifdef PRE_STEP
      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
      // FIRST (PREDICTOR) STEP
      argl.clear();
      state.set_time(time);
      state_old.set_time(time_old);
      argl.arg_add(&state,IN_VECTOR|USE_TIME_DATA);
      argl.arg_add(&state_old,IN_VECTOR|USE_TIME_DATA);
      argl.arg_add(&res,OUT_VECTOR);
      argl.arg_add(A_mom,OUT_MATRIX|PFMAT);
      argl.arg_add(&glob_param,USER_DATA);
      argl.arg_add(&wall_data,USER_DATA);
      debug.trace("-PREDICTOR- Before residual computation...");
      ierr = assemble(mesh,argl,dofmap,"comp_res_mom",&time_star);
      CHKERRA(ierr);
      debug.trace("-PREDICTOR- After residual computation.");

      do_stop = print_linear_system_and_stop && stop_on_step==tstep;
      if (do_stop && stop_mom) {
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,
				    "system.dat",&matlab); CHKERRA(ierr);
	ierr = PetscViewerSetFormat_WRAPPER(matlab,
					    PETSC_VIEWER_ASCII_MATLAB,"amom"); 
	CHKERRA(ierr);
	ierr = A_mom->view(matlab);
	ierr = PetscViewerSetFormat_WRAPPER(matlab,
					    PETSC_VIEWER_ASCII_MATLAB,"res"); 
	CHKERRA(ierr);
	ierr = VecView(res,matlab);
	ierr = PetscViewerSetFormat_WRAPPER(matlab,
					    PETSC_VIEWER_ASCII_MATLAB,"dx"); 
	CHKERRA(ierr);
	ierr = VecView(dx,matlab);
	if (!solve_system) {
	  PetscFinalize();
	  exit(0);
	}
      }

      debug.trace("-PREDICTOR- Before solving linear system...");
      ierr = A_mom->solve(res,dx); CHKERRA(ierr); 
      // Frees memory 
      ierr = A_mom->clean_mat(); CHKERRA(ierr); 
      debug.trace("-PREDICTOR- After solving linear system.");

      scal= 1.0;
      ierr = VecAXPY(&scal,dx,x);

      if (do_stop && stop_mom) {
	ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD); CHKERRA(ierr);
	PetscFinalize();
	exit(0);
      }

#if 0
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,
			     "system.dat",&matlab); CHKERRA(ierr);
      ierr = PetscViewerSetFormat_WRAPPER(matlab,
			     PETSC_VIEWER_ASCII_MATLAB,"x"); CHKERRA(ierr);
      ierr = VecView(x,matlab);
      PetscFinalize();
      exit(0);
#endif

#endif

#ifdef POI_STEP
      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
      // SECOND STEP POISSON
      ierr = VecCopy(x,xp);
      scal=0;
      ierr = VecSet(&scal,res); CHKERRA(ierr);

      if (tstep==1) {
	argl.clear();
	argl.arg_add(A_poi,OUT_MATRIX|PFMAT);
	argl.arg_add(&glob_param,USER_DATA);
	debug.trace("-POISSON- Before matrix computation...");
	ierr = assemble(mesh,argl,dofmap,"comp_mat_poi",&time_star);
	CHKERRA(ierr);
	debug.trace("-POISSON- After matrix computation.");
      }

      argl.clear();
      statep.set_time(time);	// fixme:= what time?
      argl.arg_add(&statep,IN_VECTOR|USE_TIME_DATA);
      argl.arg_add(&state,IN_VECTOR|USE_TIME_DATA);
      argl.arg_add(&res,OUT_VECTOR);
      argl.arg_add(&glob_param,USER_DATA);
      debug.trace("-POISSON- Before residual computation...");
      ierr = assemble(mesh,argl,dofmap,"comp_res_poi",&time_star);
      CHKERRA(ierr);
      debug.trace("-POISSON- After residual computation.");

#define POI_SOLVE
#ifdef POI_SOLVE
      debug.trace("-POISSON- Before solving linear system...");
      // iisdmat_stat.reset();
      ierr = A_poi->solve(res,dx); CHKERRA(ierr); 
      // iisdmat_stat.report();
      debug.trace("-POISSON- After solving linear system.");

      scal= 1.0;
      ierr = VecAXPY(&scal,dx,x);
#endif

      if (do_stop && stop_poi) {
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,
				    "system.dat",&matlab); CHKERRA(ierr);
	ierr = PetscViewerSetFormat_WRAPPER(matlab,
					    PETSC_VIEWER_ASCII_MATLAB,"apoi"); 
	CHKERRA(ierr);
	ierr = A_poi->view(matlab);
	ierr = PetscViewerSetFormat_WRAPPER(matlab,
					    PETSC_VIEWER_ASCII_MATLAB,"res"); 
	CHKERRA(ierr);
	ierr = VecView(res,matlab);
	ierr = PetscViewerSetFormat_WRAPPER(matlab,
					    PETSC_VIEWER_ASCII_MATLAB,"dx"); 
	CHKERRA(ierr);
	ierr = VecView(dx,matlab);
	PetscFinalize();
	exit(0);
      }

#endif

#ifdef PRJ_STEP
      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
      // THIRD STEP PROJECTION
      ierr = VecCopy(x,xp);
      scal=0;
      ierr = VecSet(&scal,res); CHKERRA(ierr);
      if (!reuse_mat || tstep==1) {
	argl.clear();
	statep.set_time(time);	// fixme:= what time?
	argl.arg_add(A_prj,OUT_MATRIX|PFMAT);
	debug.trace("-PROJECTION- Before matrix computation...");
	ierr = assemble(mesh,argl,dofmap,"comp_mat_prj",&time_star);
	CHKERRA(ierr);
	debug.trace("-PROJECTION- After matrix computation.");
       }
#if 0
      ierr = A_prj->view(PETSC_VIEWER_STDOUT_SELF);
#endif

      argl.clear();
      statep.set_time(time);	// fixme:= what time?
      argl.arg_add(&statep,IN_VECTOR|USE_TIME_DATA);
      argl.arg_add(&state,IN_VECTOR|USE_TIME_DATA);
      argl.arg_add(&res,OUT_VECTOR);
      if (!reuse_mat) argl.arg_add(A_prj,OUT_MATRIX|PFMAT);
      argl.arg_add(&glob_param,USER_DATA);
      debug.trace("-PROJECTION- Before residual computation...");
      ierr = assemble(mesh,argl,dofmap,"comp_res_prj",&time_star);
      CHKERRA(ierr);
      debug.trace("-PROJECTION- After residual computation.");
      if (do_stop && stop_prj) {
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,
				    "system.dat",&matlab); 
	CHKERRA(ierr);
	ierr = PetscViewerSetFormat_WRAPPER(matlab,
					    PETSC_VIEWER_ASCII_MATLAB,"aprj"); 
	CHKERRA(ierr);
      }
      
      debug.trace("-PROJECTION- Before solving linear system...");
      ierr = A_prj->solve(res,dx); CHKERRA(ierr); 
      if (!reuse_mat) { ierr = A_prj->clean_mat(); CHKERRA(ierr); }
      debug.trace("-PROJECTION- After solving linear system.");

      if (do_stop && stop_prj) {
	ierr = A_prj->view(matlab);
	ierr = PetscViewerSetFormat_WRAPPER(matlab,
					    PETSC_VIEWER_ASCII_MATLAB,"res"); 
	CHKERRA(ierr);
	ierr = VecView(res,matlab);
	ierr = PetscViewerSetFormat_WRAPPER(matlab,
					    PETSC_VIEWER_ASCII_MATLAB,"dx"); 
	CHKERRA(ierr);
	ierr = VecView(dx,matlab);
	PetscFinalize();
	exit(0);
      }

      scal= 1.0;
      ierr = VecAXPY(&scal,dx,x);
#endif
    }

    // error difference
    scal = -1.0;
    ierr = VecAXPY(&scal,x,dx_step);
    ierr  = VecNorm(dx_step,NORM_2,&norm); CHKERRA(ierr);
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
      ierr = assemble(mesh,arglf,dofmap,"gather",&time_star);
      CHKERRA(ierr);
    }

    hook_list.time_step_post(time_star.time(),tstep,gather_values);

    if (ngather>0) {
      // Print gathered values
      if (MY_RANK==0) {
	if (gather_file == "") {
	  printf("Gather results: \n");
	  for (int j=0; j < gather_values.size(); j++) 
	    printf("v_component_%d = %12.10e\n",j,gather_values[j]);
	} else {
	  gather_file_f = fopen(gather_file.c_str(),"a");
	  for (int j=0; j<gather_values.size(); j++) 
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

  if (!fractional_step) {
    delete A_tet;
    delete A_tet_c;
  } else {
    delete A_mom;
    delete A_poi;
    delete A_prj;
    ierr = VecDestroy(xp); CHKERRA(ierr); 
  }

  DELETE_SCLR(dofmap);
  DELETE_SCLR(mesh);
#ifdef DEBUG_MALLOC_USE
  fclose(malloc_log);
#endif
  PetscFinalize();
  exit(0);
}
