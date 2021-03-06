//__INSERT_LICENSE__
// $Id: absolay.cpp,v 1.2.20.1 2007/02/23 19:18:07 rodrigop Exp $
#include "./absolay.h"

AbsorbingLayer *absorbing_layer_elemset_p = NULL;
#define STEP_DBG 500
extern int abso_inwt;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AbsorbingLayer::init_hook() { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void AbsorbingLayer::initialize() {
  int ierr;
  nstep_histo = 4;

  NSGETOPTDEF_ND(int,ndim,0); //nd
  PETSCFEM_ASSERT0(ndim>0,"ndim is required");  

  // The case `ndimel=0' is a special case for concentrated
  // absorbing boundary condition, NOT distributed. 
  NSGETOPTDEF_ND(int,ndimel,-1); //nd
  if (ndimel<0) ndimel = ndim;

  NSGETOPTDEF(int,use_layer,0); //nd

  //o Include higher order term
  NSGETOPTDEF_ND(int,use_h1_term,1);

  //o If defined means replace the computation 
  //  of W for the absorbing layer by an addhoc 
  //  formula W=-ky/omega*U, with magic_abso_coef=ky/omega
  NSGETOPTDEF_ND(double,magic_abso_coef,NAN);
  NSGETOPTDEF_ND(double,magic_abso_coef1,NAN);
  PETSCFEM_ASSERT0(isnan(magic_abso_coef1)==isnan(magic_abso_coef),
                   "the magic coefs must be both defined or undefined");  
  if (!MY_RANK && !isnan(magic_abso_coef)) 
    printf("using magic coefs 0/1 %f %f\n",magic_abso_coef,magic_abso_coef1);

  NSGETOPTDEF_ND(int,nnod,-1);
  PETSCFEM_ASSERT0(nnod>0,"nnod is required");

  elem_params(nel,ndof,nelprops);
  nen = nel*ndof;

  uhist.a_resize(3,nstep_histo,nnod,ndof);
  uhist.set(0.0);
  whist.a_resize(3,nstep_histo,nnod,ndof);
  whist.set(0.0);
  u.a_resize(2,nnod,ndof);
  w.a_resize(2,nnod,ndof);

  //o Contains the Habso matrix: Hp (ndof x ndof),
  //  Hm (ndof x ndof),Uref (ndof)
  NGETOPTDEF_S(string,habso_file,"<none>");
  PETSCFEM_ASSERT0(habso_file!="<none>","habso_mat is required");  

  //o Scales volume damping term
  NSGETOPTDEF_ND(double,Kabso,1.0);
  PETSCFEM_ASSERT0(Kabso>=0.0,
                   "Kabso must be non-negatvive");  

  //o Scales volume damping term
  NSGETOPTDEF_ND(double,h1fac,NAN);
  PETSCFEM_ASSERT0(!isnan(h1fac),"h1fac is required");  

  //o Use Habso = Hm, currently only for 0-dim
  //  absorption add-hoc for examples. The absorbing matrix at
  //  the surface is computed with an add-hoc formula
  //  that is perfectly absorbing only for a given incident direction. 
  //  It's of little use in the general case, it's used here 
  //  only for validating the general absorbing layer. 
  NSGETOPTDEF_ND(int,use_addhoc_surface_abso,0);

  // If true, then the surface absorbing condition uses
  // a matrix computed for a specific incidence angle.
  // Then it is of no practical use, it is only for debugging.
  NSGETOPTDEF_ND(int,addhoc_surface_abso_use_hm,0);

  //o Characteristic time for relaxing the relation 
  //  between W and U (in secs)
  NSGETOPTDEF_ND(double,taurelax,NAN);
  PETSCFEM_ASSERT0(!isnan(taurelax),"taurelax is required");  
  
  Habso.resize(2,ndof,ndof);
  C.resize(2,ndof,ndof);
  Ay.resize(2,ndof,ndof);
  Hm.resize(2,ndof,ndof);
  Hm0.resize(2,ndof,ndof);
  Uref.resize(1,ndof);
  dvector<double> habso_data;
  habso_data.cat(habso_file.c_str());
  PETSCFEM_ASSERT(habso_data.size()==5*ndof*ndof+ndof,
                  "habso_data must be 3*ndof*ndof. ndof %d, "
                  "habso_data size %d",ndof,habso_data.size());

  Habso.set(habso_data.buff());
  C.set(habso_data.buff()+ndof*ndof);
  Ay.set(habso_data.buff()+2*ndof*ndof);
  Hm.set(habso_data.buff()+3*ndof*ndof);
  Hm0.set(habso_data.buff()+4*ndof*ndof);
  Uref.set(habso_data.buff()+5*ndof*ndof);

  H0.set(C);
  H1.set(Ay);

  if (use_addhoc_surface_abso) {
    if (addhoc_surface_abso_use_hm) Habso.set(Hm);
    else Habso.set(Hm0);
    Habso.scale(-1.0);
  }

  if (0 && !isnan(magic_abso_coef)) {
    H0.axpy(Ay,-magic_abso_coef);
    PETSCFEM_ASSERT0(!use_h1_term,"if magic_abso_coef is entered, "
                     "then use_h1_term should be deactivated");  
  }

  if (!MY_RANK) {
    printf("addhoc_surface_abso_use_hm %d\n",
           addhoc_surface_abso_use_hm);
    if (addhoc_surface_abso_use_hm) {
      printf("Using HM (not HM0) i.e. a surface absorbing condition "
             "    computed for a specific incidence angle."
             "    Then it is of no practical use, it is only for debugging.");
    }
    Habso.print("Habso: ");
    C.print("C: ");
    Ay.print("Ay: ");
    Hm.print("Hm: ");
    Hm0.print("Hm0: ");
    Uref.print("Uref:");
  }

  if (!use_addhoc_surface_abso) {
    PETSCFEM_ASSERT0(!absorbing_layer_elemset_p
                     ,"Only one absorbing layer elemset allowed");  
    absorbing_layer_elemset_p = this;
  } 

  //o Frequency for saving states
  NSGETOPTDEF(int,nsaverot,0);

  //o Frequency for saving w states
  NSGETOPTDEF_ND(int,nsaverotw,-1);
  nsaverotw = (nsaverotw>=0 ? nsaverotw : nsaverot);

  if (use_layer) {
    //o Number of elements in x direction
    NSGETOPTDEF_ND(int,Nx,-1); //nd
    PETSCFEM_ASSERT(Nx>0,"Nx is required. Nx %d",Nx);  
    
    //o Number of elements in y direction
    NSGETOPTDEF_ND(int,Ny,-1); //nd
    PETSCFEM_ASSERT(Ny>0,"Ny is required. Ny %d",Ny);  
    
    //o Size of elements in y direction
    NSGETOPTDEF_ND(double,hy,NAN); //nd
    PETSCFEM_ASSERT(!isnan(hy),"hy is required. hy %f",hy);  
    PETSCFEM_ASSERT(hy>0.0,"hy must be positive. hy %f",hy);  

    //o Size of elements in x direction
    NSGETOPTDEF_ND(double,hx,NAN); //nd
    PETSCFEM_ASSERT(!isnan(hx),"hx is required. hx %f",hx);  
    PETSCFEM_ASSERT(hx>0.0,"hx must be positive. hx %f",hx);  
  }

  //o Time step
  NSGETOPTDEF_ND(double,Dt,NAN); //nd
  PETSCFEM_ASSERT(!isnan(Dt),"Dt is required. Dt %f",Dt);  
  PETSCFEM_ASSERT(Dt>0.0,"Dt must be positive. Dt %f",Dt);  

  kfac = 1.0/(1.0+2.0*Dt/taurelax/3.0);
  // kfac = 1.0; // equivalent to taurelax = Inf
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "int advective::ask(char *,int &)"
int AbsorbingLayer::ask(const char *jobinfo,int &skip_elemset) {
   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_prof);
   return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
double rnd(double) {
  return 2.0*drand()-1.0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void AbsorbingLayer::node2jk(int node,int &j,int &k) {
  j = node/(Ny+1);
  k = node%(Ny+1);
  k = modulo(k,Ny);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int AbsorbingLayer::jk2node(int j,int k) {
  int kk = k;
  if (kk<0 || kk>=Ny+1) kk = modulo(kk,Ny+1);
  return (Ny+1)*j+kk;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void AbsorbingLayer::
new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
	     const Dofmap *dofmap,const char *jobinfo,
	     const ElementList &elemlist,
	     const TimeData *time_data) try {

  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_prof);

  int ierr=0;

  int locdof,kldof,lldof;

  // Unpack Dofmap
  if (nnod!=dofmap->nnod) {
    printf("nnod %d != dofmap->nnod %d\n",nnod,dofmap->nnod);
    abort();
  }

  // Unpack nodedata
  int nu=nodedata->nu;
  // H is a generalized local property passed per node with the nodal
  // coordinates. In shallow water nH =1 and H is the depth. It is
  // needed in order to compute the source term. In 1D Euler it may be
  // the area section of the tube. Its gradient is needed for the
  // source term in the momentum eqs.
  int nH = nu-ndim;
  FMatrix  Hloc(nel,nH),H(nH);
  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  // lambda_max:= the maximum eigenvalue of the jacobians.
  // used to compute the critical time step.
  vector<double> *dtmin;
  int PFUNUSED jdtmin;
  GlobParam *glob_param=NULL;
  double alpha = 1.0, alpha_glob=NAN;
  // The trapezoidal rule integration parameter
  arg_data PFUNUSED *staten=NULL,*stateo=NULL,*retval=NULL,
    *fdj_jac=NULL,*jac_prof=NULL,*Ajac=NULL;
  if (comp_res) {
    int j=-1;
    stateo = &arg_data_v[++j]; //[0]
    staten = &arg_data_v[++j]; //[1]
    retval  = &arg_data_v[++j];//[2]
    jdtmin = ++j;//[3]
#define DTMIN ((*(arg_data_v[jdtmin].vector_assoc))[0])
#define WAS_SET arg_data_v[jdtmin].was_set
    Ajac = &arg_data_v[++j];//[4]
    glob_param = (GlobParam *)arg_data_v[++j].user_data;;
    alpha_glob = (glob_param->alpha);

    if (ADVDIF_CHECK_JAC)
      fdj_jac = &arg_data_v[++j];
  }

  FastMat2 matloc(4,nel,ndof,nel,ndof);

  //o Uses operations caches for computations with the FastMat2
  //  library for internal computations. Note that this affects also the
  //  use of caches in routines like fluxes, etc...
  NSGETOPTDEF(int,use_fastmat2_cache,1);

  // It seems that when using penalization we must 
  // use alpha=1.0, because otherwise it is unstable. 
  // However due to how the Newton scheme
  // is implemented in advdif, we must return in matloc 
  // the Jacobian as matloc = -d(veccontr)/dU/alpha_glob
  // where alpha_glob is the alpha used in the program. 
  alpha = 1.0;
  if (!use_layer) ndimel = 0;

  matloc.set(1.0);
  if (comp_prof) {
    jac_prof = &arg_data_v[0];
    matloc.export_vals(jac_prof->profile);
  }

  int PFUNUSED nlog_vars;
  const int *log_vars;
  // adv_diff_ff->get_log_vars(nlog_vars,log_vars);
  //o Use log-vars for $k$ and $\epsilon$
  NSGETOPTDEF(int,use_log_vars,0);
  if (!use_log_vars) nlog_vars=0;

  // Allocate local vecs
  FMatrix veccontr(nel,ndof),veccontr_mass(nel,ndof),
    xloc(nel,ndim),lstate(nel,ndof),
    lstateo(nel,ndof),lstaten(nel,ndof),dUloc_c(nel,ndof),
    dUloc(nel,ndof);
  FastMat2 Un(1,ndof),Uo(1,ndof),Ualpha(1,ndof),
    Jaco(2,ndimel,ndim),iJaco(2,ndimel,ndimel),
    dU(1,ndof),dshapex(2,ndimel,nel), normal(1,ndim),tmp1,tmp2,
    shape(1,nel),tmp3,W(1,ndof),Wbar(1,ndof),W1(1,ndof),W2(1,ndof),
    Wmagic(1,ndof),U1(1,ndof),dudy(1,ndof);

  nen = nel*ndof;

  //o Type of element geometry to define Gauss Point data
  NGETOPTDEF_S(string,geometry,cartesian2d);
  GPdata gp_data;
  int npg=1;
  if (0 && use_layer>0) {
    NSGETOPTDEF_ND(int,npg,0); //nd
    assert(npg>0);
    gp_data.init(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);
  }

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE	 (*gp_data.FM2_shape[ipg])
#define WPG	 (gp_data.wpg[ipg])

  double detJaco, wpgdet;
  int elem, ipg,node, jdim, kloc,lloc,ldof;

  // Position of current element in elemset and in chunk
  int k_elem, k_chunk;

  FastMatCacheList cache_list;
  if (use_fastmat2_cache) FastMat2::activate_cache(&cache_list);
  vector<int> locnodes(nel);
  double kvol = Kabso*hy*hx;

  for (ElementIterator element = elemlist.begin();
       element!=elemlist.end(); element++) try {

    element.position(k_elem,k_chunk);
    FastMat2::reset_cache();

    // Initialize element
    // adv_diff_ff->element_hook(element);

    // Get nodedata info (coords. etc...)
    element.node_data(nodedata,xloc.storage_begin(),
		      Hloc.storage_begin());

    if (comp_prof) {
      matloc.export_vals(element.ret_mat_values(*jac_prof));
      continue;
    }

    if (comp_res) {
      lstateo.set(element.vector_values(*stateo));
      lstaten.set(element.vector_values(*staten));
    }

    // State at time t_{n+\alpha}
    lstate.set(0.).axpy(lstaten,alpha).axpy(lstateo,(1-alpha));
    veccontr.set(0.);
    matloc.set(0.0);

    for (ipg=0; ipg<npg; ipg++) {
      
      if (0 && use_layer>0) {
        //      Matrix xpg = SHAPE * xloc;
        Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);
        
        PETSCFEM_ASSERT0(ndim==ndimel,"not implemented yet");
        // iJaco.inv(Jaco);
        detJaco = Jaco.det();
        
        if (detJaco<=0.) {
          int k,ielh;
          element.position(k,ielh);
          detj_error(detJaco,k);
          set_error(1);
        }
        wpgdet = detJaco*WPG;
        shape.set(SHAPE);
      } else {
        shape.set(1.0);
        wpgdet = 1.0;
      }

      assert(!isnan(alpha_glob));
      if (nel==1) {
        Un.prod(shape,lstaten,-1,-1,1);
        Uo.prod(shape,lstateo,-1,-1,1);
        Ualpha.set(0.).axpy(Uo,1-alpha).axpy(Un,alpha);
        tmp1.prod(Habso,shape,shape,2,4,1,3);
        matloc.axpy(tmp1,alpha*Kabso*wpgdet/alpha_glob);
        dU.set(Ualpha).minus(Uref);
        tmp2.prod(shape,Habso,dU,1,2,-1,-1);
        veccontr.axpy(tmp2,-Kabso*wpgdet);
      } else if (nel==3) {
        element_connect(element,&locnodes[0]);
        int node = locnodes[1];
        // int j,k,lN,lS;
        // node2jk(node-1,j,k);
        // lN = jk2node(j,k+1);
        // lS = jk2node(j,k-1);
        // assert(lN==locnodes[2]-1);
        // assert(lS==locnodes[0]-1);
        W1.set(&whist.e(0,node-1,0));
        W2.set(&whist.e(1,node-1,0));
        U1.set(&uhist.e(0,node-1,0));
        Wbar.set(0.0).axpy(W1,4.0).axpy(W2,-1.0);

        veccontr.set(0.0);
        lstate.ir(1,2);
        dU.set(lstate).minus(Uref);
        veccontr.ir(1,2).prod(H0,dU,1,-1,-1);
        int use_wmagic = 0;
        if (use_h1_term) {
          lstate.ir(1,3);
          W.set(lstate);
          lstate.ir(1,1);
          W.minus(lstate).scale(0.5/hy);
          dudy.set(W);
          W.scale(2.0*Dt).add(Wbar).scale(kfac/3.0);
          // if (0 && abso_current_step==STEP_DBG && abso_inwt==1)
          lstate.rs();
          lstate.ir(1,2);
          Wmagic.set(lstate).scale(magic_abso_coef)
            .axpy(U1,-magic_abso_coef1);
          lstate.rs();
#define NODE_DBG 1662
          if (0 && node==NODE_DBG && abso_inwt==1) {
            double 
              *Up = dU.storage_begin(),
              *Wp = W.storage_begin(),
              *Wmagicp = Wmagic.storage_begin(),
              *Wbarp = Wbar.storage_begin(),
              *dudyp = dudy.storage_begin();
            printf("node %d step %d U %f %f W %f %f Wmagic "
                   "%f %f dudy %f %f Wbar %f %f\n",
                   NODE_DBG,abso_current_step,
                   Up[0],Up[1],Wp[0],Wp[1],Wmagicp[0],Wmagicp[1],
                   dudyp[0],dudyp[1],Wbarp[0],Wbarp[1]);
          }
          if (!isnan(magic_abso_coef)) {
            if (print_w_values && abso_inwt==1) {
              lstate.rs();
              double *p;
              nodes.push_back(node);
              p = dU.storage_begin();
              values.push_back(p[0]);
              values.push_back(p[1]);
              p = W.storage_begin();
              values.push_back(p[0]);
              values.push_back(p[1]);
              p = Wmagic.storage_begin();
              values.push_back(p[0]);
              values.push_back(p[1]);
            }
          }
          if (use_wmagic) tmp3.prod(H1,Wmagic,1,-1,-1);
          else tmp3.prod(H1,W,1,-1,-1);
          tmp3_max = max(tmp3.norm_p_all(),tmp3_max);
          veccontr.axpy(tmp3,h1fac).rs();
        }

        matloc.rs().set(0.0).ir(1,2);
        matloc.ir(3,2).axpy(H0,1.0);
        if (use_h1_term) {
          if (use_wmagic) {
            matloc.axpy(Ay,-magic_abso_coef);
          } else {
            double cfac = h1fac*kfac*Dt/(3.0*hy);
            matloc.ir(3,1).axpy(H1,-cfac);
            matloc.ir(3,3).axpy(H1,+cfac);
          }
        }

        veccontr.rs().scale(-kvol);
        matloc.rs().scale(kvol/alpha_glob);
      }
    }

    veccontr.export_vals(element.ret_vector_values(*retval));
    matloc.export_vals(element.ret_mat_values(*Ajac));
    
  } catch (GenericError e) {
    set_error(1);
    return;
  }

  FastMat2::void_cache();
  FastMat2::deactivate_cache();
} catch (GenericError e) {
  set_error(1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AbsorbingLayer::time_step_pre(int step) { 
  tmp3_max = 0.0;
  if (!MY_RANK && step==STEP_DBG) 
    printf("in pre: entering step %d\n",step);
  abso_current_step=step;
  print_w_values = (step%50==0);
  print_w_values = 0;
  if (print_w_values) {
    nodes.clear();
    values.clear();
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AbsorbingLayer::time_step_post(int step) {
  double tmp;
  MPI_Reduce(&tmp3_max,&tmp,1,MPI_DOUBLE,MPI_MAX,
             0,PETSCFEM_COMM_WORLD); 
  tmp3_max = tmp;
  if (0 && !MY_RANK) {
    printf("in abso_hook::time_step_post()\n");
    printf("max ||tmp3|| %f\n",tmp3_max);
  }
  u.read("./fsabso2d.state.tmp");

  for (int j=nstep_histo-1; j>=1; j--) {
    for (int l=0; l<nnod; l++) {
      for (int k=0; k<ndof; k++) {
        uhist.e(j,l,k) = uhist.e(j-1,l,k);
        whist.e(j,l,k) = whist.e(j-1,l,k);
      }
    }
  }

  for (int l=0; l<nnod; l++) 
    for (int k=0; k<ndof; k++) 
      uhist.e(0,l,k) = u.e(l,k);

  double DUDY[2]={NAN,NAN};
  for (int l=0; l<nnod; l++) {
    // int 
    // j = l/(Ny+1),             // x-position of node l
    // kk = l%(Ny+1),            // y-position of node l
    // kN,                       // y-position of North node
    // kS,                       // y-position of North node
    // lN,                       // node at North of l
    // lS;                       // node at North of l
    // kN = kk+1;
    // if (kN>=Ny+1) kN -= Ny;
    // kS = kk-1;
    // if (kS<0) kS += Ny;
    // lN = (Ny+1)*j+kN;           // node at North of l
    // lS = (Ny+1)*j+kS;           // node at North of l

    int j,k,lN,lS,kk;
    node2jk(l,j,k);
    lN = jk2node(j,modulo(k+1,Ny));
    lS = jk2node(j,modulo(k-1,Ny));
    if (0 && !MY_RANK && l==NODE_DBG-1) {
      printf("lS %d, l %d, lN %d\n",lS,l,lN);
      printf("uS %f uN %f\n",uhist.e(0,lS,0),uhist.e(0,lN,0));
    }

    for (int kdof=0; kdof<ndof; kdof++) {
      double 
        dudy = (uhist.e(0,lN,kdof)-uhist.e(0,lS,kdof))/(2.0*hy),
        ww = kfac*(2.0*Dt*dudy+4.0*whist.e(1,l,kdof)
                   -whist.e(2,l,kdof))/3.0;
      if (l==NODE_DBG-1) DUDY[kdof] = dudy;
      whist.e(0,l,kdof) = ww;
      w.e(l,kdof) = ww;
    }
  }
  if (0 && !MY_RANK) {
    printf("node %d dudy %f %f uS %f uO %f uN %f uprev %f wh %f %f\n",
           NODE_DBG,DUDY[0],DUDY[1],
           uhist.e(0,NODE_DBG-2,0),
           uhist.e(0,NODE_DBG-1,0),
           uhist.e(0,NODE_DBG-0,0),
           uhist.e(1,NODE_DBG-1,0),
           whist.e(0,NODE_DBG-1,0),
           whist.e(1,NODE_DBG-1,0));
    printf("node %d W0 %f %f W1 %f %f W2 %f %f\n",
           NODE_DBG,
           whist.e(0,NODE_DBG-1,0),
           whist.e(0,NODE_DBG-1,1),
           whist.e(1,NODE_DBG-1,0),
           whist.e(1,NODE_DBG-1,1),
           whist.e(2,NODE_DBG-1,0),
           whist.e(2,NODE_DBG-1,1));
  }

  if (0 && !MY_RANK && step==STEP_DBG) {
    w.print("w.state.tmp");
    u.print("u.state.tmp");
  }

  // if (0 && step==(STEP_DBG+2)) {
  if (print_w_values) {
    if(!MY_RANK) 
      PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,
                              "-------- step %d: node u w wmagic ----------\n",
                              abso_current_step);
    int ndata=6;
    int OKG,OK = (values.size()==ndata*nodes.size());
    MPI_Allreduce(&OK,&OKG,1,MPI_INT,MPI_LAND,PETSCFEM_COMM_WORLD);
    if (!OKG) {
      PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,
                              "[%d] values %d nodes %d ndata %d OK %d\n",
                              MY_RANK,values.size(),nodes.size(),ndata,OK);
      PetscSynchronizedFlush(PETSCFEM_COMM_WORLD); 
      PETSCFEM_ERROR0("Not OK");  
    }
    for (unsigned int j=0; j<nodes.size(); j++) {
      PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,"%d ",nodes[j]);
      for (int k=0; k<ndata; k++) 
        PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,"%f ",values[j*ndata+k]);
      PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,"\n");
      
    }
    PetscSynchronizedFlush(PETSCFEM_COMM_WORLD); 
  }

  if (!MY_RANK && nsaverotw>0 && (step % nsaverotw == 0)) {
    char file[100];
    int indx = step/nsaverotw-1;
    sprintf(file,"./STEPS/fsabso2d.w-%d.tmp",indx);
    w.print(file);
  }
}
