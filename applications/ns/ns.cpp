//__INSERT_LICENSE__
//$Id: ns.cpp,v 1.73 2002/04/04 01:42:17 mstorti Exp $

//#define ROCKET_MODULE 
#ifndef ROCKET_MODULE 
#include <src/debug.h>
#include <malloc.h>

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/sttfilter.h>
#include <src/pfmat.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "elast.h"
#include "qharm.h"
#include "gatherer.h"

#include <applications/ns/nsi_tet.h>
#include <applications/ns/nssup.h>

static char help[] = "PETSc-FEM Navier Stokes module\n\n";

extern int MY_RANK,SIZE;
TextHashTable *GLOBAL_OPTIONS;

//debug:=
int TSTEP=0;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "bless_elemset"
void bless_elemset(char *type,Elemset *& elemset) {
  //  SET_ELEMSET_TYPE(internal)
  //    SET_ELEMSET_TYPE(fracstep) Por ahora lo desactivamos hasta que
  // hagamos la interfase
  // SET_ELEMSET_TYPE(nsi_tet)
  //  SET_ELEMSET_TYPE(nsi_tet_les)
    SET_ELEMSET_TYPE(qharm)
    SET_ELEMSET_TYPE(ns_id)
    SET_ELEMSET_TYPE(ns_sup)
    SET_ELEMSET_TYPE(ns_sup_res)
      // SET_ELEMSET_TYPE(elasticity_f)
    SET_ELEMSET_TYPE(elasticity)
    SET_ELEMSET_TYPE(nsi_tet_les_fm2)
    SET_ELEMSET_TYPE(nsi_tet_les_ther)
    SET_ELEMSET_TYPE(nsi_tet_keps)
    SET_ELEMSET_TYPE(nsi_tet_keps_rot)
    SET_ELEMSET_TYPE(bcconv_ns_fm2)
    SET_ELEMSET_TYPE(bcconv_nsther_fm2)
    SET_ELEMSET_TYPE(wall)
    SET_ELEMSET_TYPE(wallke)
    SET_ELEMSET_TYPE(wall_law_res)
    SET_ELEMSET_TYPE(force_integrator)
	{
	printf("not known elemset \"type\": %s\n",type);
	exit(1);
	}
}


//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
#undef __FUNC__
#define __FUNC__ "main"
int main(int argc,char **args) {

  Vec     x, dx, xold,
    dx_step, res;		// approx solution, RHS, residual
  Viewer matlab;
  PFMat *A_tet, *A_tet_c;			// linear system matrix 
  double  norm, *sol, scal;	// norm of solution error
  int     ierr, i, n = 10, col[3], flg, size, node,
    jdof, k, kk, nfixa,
    kdof, ldof, lloc, nel, nen, neq, nu,
    myrank;
  // Initialize time
  Time time,time_old,time_star; time.set(0.);
  GlobParam glob_param;
  string save_file_res;
  State state(x,time),state_old(xold,time_old);

  // ierr = MatCreateShell(PETSC_COMM_WORLD,int m,int n,int M,int N,void *ctx,Mat *A)
  char fcase[FLEN+1];
  Dofmap *dofmap = new Dofmap;
  Mesh *mesh;
  // arglf:= argument list for computing gathered quantities as forces
  arg_list argl, arglf;
  vector<double> hmin;
#ifdef RH60   // fixme:= STL vector compiler bug??? see notes.txt
  hmin.resize(1);
#endif

  PetscInitialize(&argc,&args,(char *)0,help);
  print_copyright();
  PetscPrintf(PETSC_COMM_WORLD,"-------- Navier-Stokes module ---------\n");

  // Get MPI info
  MPI_Comm_size(PETSC_COMM_WORLD,&SIZE);
  MPI_Comm_rank(PETSC_COMM_WORLD,&MY_RANK);

  //  if (size != 1) SETERRA(1,0,"This is a uniprocessor example only!");
  ierr = OptionsGetString(PETSC_NULL,"-case",fcase,FLEN,&flg);
  CHKERRA(ierr);
  if (!flg) {
    PetscPrintf(PETSC_COMM_WORLD,
		"Option \"-case <filename>\" not passed to PETSc-FEM!!\n");
    PetscFinalize();
    exit(0);
  }

  // Read the mesh
  read_mesh(mesh,fcase,dofmap,neq,SIZE,MY_RANK);

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

  //o Dimension of the problem.
  GETOPTDEF(int,ndim,3);

  PetscPrintf(PETSC_COMM_WORLD,"After readmesh...\n");

  //o Number of inner iterations for the global non-linear
  // Newton  problem. 
  GETOPTDEF(int,nnwt,1);
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
  //o Update jacobian each $n$-th time step. 
  GETOPTDEF(int,update_jacobian_steps,1);
  assert(update_jacobian_steps>=1);
  //o Update jacobian each $n$-th time step. 
  GETOPTDEF(int,update_jacobian_start_steps,INF);
  assert(update_jacobian_start_steps>=0);
#undef INF

  //o _T: vector<int>
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
  
  GETOPTDEF(int,verify_jacobian_with_numerical_one,0);
  //o After computing the linear system solves it and prints Jacobian,
  // right hand side and solution vector, and stops. 
  GETOPTDEF(int,print_linear_system_and_stop,0);
  //o Print the residual each \verb+nsave+ steps. 
  GETOPTDEF(int,print_residual,0);
  //o Solve system before \verb+print\_linear_system_and_stop+
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
  //o Flag if steady solution or not (uses Dt=inf). If \verb+steady+
  // is set to 1, then the computations are as if $\Dt=\infty$. 
  // The value of \verb+Dt+ is used for printing etc... If \verb+Dt+
  // is not set and \verb+steady+ is set then \verb+Dt+ is set to one.
  GETOPTDEF(int,steady,0);
  if (steady && Dt==0.) Dt=1.;
  // Set values to be passed through global options
  glob_param.steady=steady;
  glob_param.Dt=Dt;
  //o Trapezoidal method parameter. \verb+alpha=1+:
  // Backward Euler. \verb+alpha=0+: Forward Euler.
  // \verb+alpha=0.5+: Crank-Nicholson. 
  GETOPTDEF(double,alpha,1.); 
  glob_param.alpha=alpha;
  //o Number of ``gathered'' quantities.
  GETOPTDEF(int,ngather,0);
  //o Print values in this file 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,gather_file,"gather.out");
  // Initialize gather_file
  FILE *gather_file_f;
  if (MY_RANK==0 && ngather>0) {
    gather_file_f = fopen(gather_file.c_str(),"w");
    fprintf(gather_file_f,"");
    fclose(gather_file_f);
  }

  //o Use the LES/Smagorinsky turbulence model. 
  GETOPTDEF(int,LES,0);

  //o Use IISD (Interface Iterative Subdomain Direct) or not.
  GETOPTDEF(int,use_iisd,0);
  //o Type of solver. May be \verb+iisd+ or \verb+petsc+. 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,solver,petsc);
  if (use_iisd) solver = string("iisd");

  //o The pattern to generate the file name to save in for
  // the rotary save mechanism.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file_pattern,outvector%d.out);

  //o The name of the file to save the state vector. 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file,outvector.out);
  save_file_res = save_file + string(".res");

  //o Name of file where to read the nodes for the ``print some'' 
  // feature. 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,print_some_file,);
  //o Name of file where to save node values for the ``print some'' 
  // feature. 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file_some,outvsome.out);
  //o Print, after execution, a report of the times a given option
  // was accessed. Useful for detecting if an option was used or not.
  GETOPTDEF(int,report_option_access,1);

  set<int> node_list;
  print_some_file_init(mesh->global_options,
		       print_some_file.c_str(),
		       save_file_some.c_str(),node_list);

  // initialize vectors
  dofmap->create_MPI_vector(x);
  State sx(x,time);		// convert to State

  LPFilterGroup filter(GLOBAL_OPTIONS,sx,Dt);
  
  // Use IISD (Interface Iterative Subdomain Direct) or not.
  // A_tet = (use_iisd ? &IISD_A_tet : &PETSc_A_tet);
  A_tet = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());
  A_tet_c = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());

#if 0
  const int NT=200;
  const double DT=0.1, omega=1.;
  for (int kk=0; kk<NT; kk++) {
    double t=double(kk)*DT;
    double a=1.+sin(omega*t);
    sx.set_cnst(a);
    ff.update();
    PetscPrintf(PETSC_COMM_WORLD,"t = %f,  a = %f\n",t,a);
    //print_some("filter.some",ff.state(),dofmap,node_list);
    print_some("filter.some.tmp",sx,dofmap,node_list);
    print_some("filter.some.tmp",ff.state(),dofmap,node_list);
  }
  PetscFinalize();
  exit(0);
#endif
  
  ierr = VecDuplicate(x,&xold); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx_step); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx); CHKERRA(ierr);
  ierr = VecDuplicate(x,&res); CHKERRA(ierr);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // initialize state vectors
  scal=0;
  ierr = VecSet(&scal,x); CHKERRA(ierr);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // COMPUTE ACTIVE PROFILE
  VOID_IT(argl);
  argl.arg_add(A_tet,PROFILE|PFMAT);
  debug.trace("Computing profile...");
  ierr = assemble(mesh,argl,dofmap,"comp_mat",&time); CHKERRA(ierr); 
  debug.trace("After computing profile.");

#if 0 //dbg
  VOID_IT(argl);
  argl.arg_add(A_tet,OUT_MATRIX);
  for (int jjj=0; jjj<10; jjj++) {
    printf("[loop iter %d]\n",jjj);
    ierr = assemble(mesh,argl,dofmap,"comp_mat",&time); CHKERRA(ierr); 
  }
  PetscFinalize();
  exit(0);
#endif

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Build octree for nearest neighbor calculation

  vector<double> *data_pts_ = new vector<double>;
  vector<ElemToPtr> *elemset_pointer = new vector<ElemToPtr>;
  WallData *wall_data;
  if (LES) {
    VOID_IT(argl);
    argl.arg_add(data_pts_,USER_DATA);
    argl.arg_add(elemset_pointer,USER_DATA);
    Elemset *elemset=NULL;
    argl.arg_add(elemset,USER_DATA);
    ierr = assemble(mesh,argl,dofmap,"build_nneighbor_tree",&time); CHKERRA(ierr); 
    PetscPrintf(PETSC_COMM_WORLD,"After nearest neighbor tree.\n");

    wall_data = new WallData(data_pts_,elemset_pointer,ndim);

    //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    // Find nearest neighbor for each volume element
    VOID_IT(argl);
    argl.arg_add(wall_data,USER_DATA);
    ierr = assemble(mesh,argl,dofmap,"get_nearest_wall_element",
		    &time); CHKERRA(ierr); 
  }

#if 0 
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  ierr = opt_read_vector(mesh,x,dofmap,MY_RANK); CHKERRA(ierr);

  // Filter *filter(x,*mesh);

  // update_jacobian_this_step:= Flags whether this step the
  // jacobian should be updated or not 
  int update_jacobian_this_step,update_jacobian_this_iter;
  for (int tstep=1; tstep<=nstep; tstep++) {
    TSTEP=tstep; //debug:=
    time_old.set(time.time());
    time_star.set(time.time()+alpha*Dt);
    time.inc(Dt);
    PetscPrintf(PETSC_COMM_WORLD,"Time step: %d, time: %g %s\n",
		tstep,time.time(),(steady ? " (steady) " : ""));
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    // TET ALGORITHM

    // Jacobian update logic
    update_jacobian_this_step = (tstep < update_jacobian_start_steps) 
      || ((tstep-update_jacobian_start_steps) % update_jacobian_steps == 0);
    
    // Inicializacion del paso
    ierr = VecCopy(x,dx_step);
    ierr = VecCopy(x,xold);

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
	// ierr = A_tet->clean_factor(); CHKERRA(ierr); 
      }

      // Compute wall stresses
      VOID_IT(argl);
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
      if (LES && SIZE>1) {
	VOID_IT(argl);
	ierr = assemble(mesh,argl,dofmap,"communicate_shear_vel",
			&time_star); CHKERRA(ierr);
      }

      scal=0;
      ierr = VecSet(&scal,res); CHKERRA(ierr);
      if (update_jacobian_this_iter) {
	ierr = A_tet->clean_mat(); CHKERRA(ierr); 
      }

      VOID_IT(argl);
      state.set_time(time);
      state_old.set_time(time_old);
      argl.arg_add(&state,IN_VECTOR|USE_TIME_DATA);
      argl.arg_add(&state_old,IN_VECTOR|USE_TIME_DATA);
      argl.arg_add(&res,OUT_VECTOR);
      if (update_jacobian_this_iter) argl.arg_add(A_tet,OUT_MATRIX|PFMAT);
#ifdef RH60 // fixme:= STL vector compiler bug??? see notes.txt
      argl.arg_add(&hmin,VECTOR_MIN);
#else
      argl.arg_add(&hmin,USER_DATA);
#endif
      argl.arg_add(&glob_param,USER_DATA);
      argl.arg_add(wall_data,USER_DATA);

      const char *jobinfo = (update_jacobian_this_iter ? "comp_mat_res" : "comp_res");

      // In order to measure performance
      if (measure_performance) {
	ierr = measure_performance_fun(mesh,argl,dofmap,jobinfo,
			&time_star); CHKERRA(ierr);
	PetscFinalize();
	exit(0);
      }

      debug.trace("Before residual computation...");
      ierr = assemble(mesh,argl,dofmap,jobinfo,&time_star); CHKERRA(ierr);
      debug.trace("After residual computation.");

#if 0
      ierr = ViewerASCIIOpen(PETSC_COMM_WORLD,
 			     "system.dat",&matlab); CHKERRA(ierr);
      ierr = ViewerSetFormat(matlab,
			     VIEWER_FORMAT_ASCII_MATLAB,"a_tet"); CHKERRA(ierr);
      ierr = A_tet->view(matlab); CHKERRA(ierr); 
#endif

#if 0 //debug:=
      ierr = ViewerASCIIOpen(PETSC_COMM_WORLD,
			     "system.dat",&matlab); CHKERRA(ierr);
      ierr = ViewerSetFormat(matlab,
			     VIEWER_FORMAT_ASCII_MATLAB,"x"); CHKERRA(ierr);
      ierr = VecView(x,matlab);
      ierr = ViewerSetFormat(matlab,
			     VIEWER_FORMAT_ASCII_MATLAB,"res"); CHKERRA(ierr);
      ierr = VecView(res,matlab);
      PetscFinalize();
      exit(0);
#endif

      Viewer matlab;
      if (verify_jacobian_with_numerical_one) {
	ierr = ViewerASCIIOpen(PETSC_COMM_WORLD,
			       "system.dat.tmp",&matlab); CHKERRA(ierr);
	ierr = ViewerSetFormat(matlab,
			       VIEWER_FORMAT_ASCII_MATLAB,
			       "atet"); CHKERRA(ierr);

	ierr = A_tet->view(matlab); CHKERRQ(ierr); 
	
	ierr = A_tet_c->duplicate(MAT_DO_NOT_COPY_VALUES,*A_tet); CHKERRA(ierr);
	ierr = A_tet->clean_mat(); CHKERRA(ierr); 
	ierr = A_tet_c->clean_mat(); CHKERRA(ierr); 

	VOID_IT(argl);
	argl.arg_add(&x,PERT_VECTOR);
	argl.arg_add(&xold,IN_VECTOR);
	argl.arg_add(A_tet_c,OUT_MATRIX_FDJ|PFMAT);
	if (update_jacobian_this_iter) argl.arg_add(A_tet,OUT_MATRIX|PFMAT);
#ifdef RH60    // fixme:= STL vector compiler bug??? see notes.txt
	argl.arg_add(&hmin,VECTOR_MIN);
#else
	argl.arg_add(&hmin,USER_DATA);
#endif
	argl.arg_add(&glob_param,USER_DATA);
	argl.arg_add(wall_data,USER_DATA);
	ierr = assemble(mesh,argl,dofmap,jobinfo,
			&time_star); CHKERRA(ierr);

	ierr = ViewerSetFormat(matlab,
			       VIEWER_FORMAT_ASCII_MATLAB,"atet_fdj"); CHKERRA(ierr);
	ierr = A_tet_c->view(matlab); CHKERRQ(ierr); 

	PetscFinalize();
	exit(0);
      }

      // A_tet->build_sles(GLOBAL_OPTIONS);

      if (!print_linear_system_and_stop || solve_system) {
	// ierr = SLESSolve(sles_tet,res,dx,&its); CHKERRA(ierr); 
	debug.trace("Before solving linear system...");
	ierr = A_tet->solve(res,dx); CHKERRA(ierr); 
	debug.trace("After solving linear system.");
      }

      if (print_linear_system_and_stop) {
	PetscPrintf(PETSC_COMM_WORLD,
		    "Printing residual and matrix for"
		    " debugging and stopping.\n");
	ierr = ViewerASCIIOpen(PETSC_COMM_WORLD,
			       "system.dat",&matlab); CHKERRA(ierr);
	ierr = ViewerSetFormat(matlab,
			       VIEWER_FORMAT_ASCII_MATLAB,"atet"); CHKERRA(ierr);
	ierr = A_tet->view(matlab);
	ierr = ViewerSetFormat(matlab,
			       VIEWER_FORMAT_ASCII_MATLAB,"res"); CHKERRA(ierr);
	ierr = VecView(res,matlab);

	if (solve_system) {
	  ierr = ViewerSetFormat(matlab,
				 VIEWER_FORMAT_ASCII_MATLAB,"dx");
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
      
      scal= relfac/alpha;
      ierr = VecAXPY(&scal,dx,x);

#if 0
      ierr = VecView(x,VIEWER_STDOUT_WORLD); CHKERRA(ierr);
      PetscFinalize();
      exit(0);
#endif

      // fixme:= SHOULD WE CHECK HERE FOR NEWTON CONVERGENCE?

    } // end of loop over Newton subiteration (inwt)

    // error difference
    scal = -1.0;
    ierr = VecAXPY(&scal,x,dx_step);
    ierr  = VecNorm(dx_step,NORM_2,&norm); CHKERRA(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"============= delta_u = %10.3e\n",norm);
    print_vector_rota(save_file_pattern.c_str(),x,dofmap,&time,
		      tstep-1,nsaverot,nrec,nfile);
  
    if (tstep % nsave == 0) {
      print_vector(save_file.c_str(),x,dofmap,&time);
      if (print_residual) 
	print_vector(save_file_res.c_str(),res,dofmap,&time);
    }

    if (normres_external < tol_newton) {
      PetscPrintf(PETSC_COMM_WORLD,
		  "Tolerance on newton loop reached:  "
		  "|| R ||_0,  norm_res =%g < tol = %g\n",
		  normres_external,tol_newton);
      break;
    }

    // Compute gathered quantities, for instance total force on walls
    if (ngather>0) {
      vector<double> gather_values(ngather);
      arglf.clear();
      arglf.arg_add(&state,IN_VECTOR|USE_TIME_DATA);
      arglf.arg_add(&state_old,IN_VECTOR|USE_TIME_DATA);
      arglf.arg_add(&gather_values,VECTOR_ADD);
      ierr = assemble(mesh,arglf,dofmap,"gather",&time_star);
      CHKERRA(ierr);
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

  if (report_option_access && MY_RANK==0) TextHashTable::print_stat();
  print_vector(save_file.c_str(),x,dofmap,&time);

  // ierr = VecDestroy(x); CHKERRA(ierr); 
  ierr = VecDestroy(dx); CHKERRA(ierr); 
  ierr = VecDestroy(res); CHKERRA(ierr); 

  delete A_tet;
  delete A_tet_c;

#ifdef DEBUG_MALLOC_USE
  fclose(malloc_log);
#endif
  PetscFinalize();
  exit(0);
}

#endif
