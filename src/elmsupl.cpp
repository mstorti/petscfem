//__INSERT_LICENSE__
//$Id: elmsupl.cpp,v 1.21 2003/09/14 13:20:21 mstorti Exp $

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <vector>
#include <set>

#include <src/fem.h>
#include <src/utils.h>
#include <src/getprop.h>
#include <src/elemset.h>
#include <src/idmap.h>
#include <src/dofmap.h>
#include <src/arglist.h>
#include <src/readmesh.h>
#include <src/pfmat.h>
#include <src/timestat.h>
#include <src/util3.h>
#include <src/autostr.h>
#include <src/dvector.h>

// iteration modes
#define NOT_INCLUDE_GHOST_ELEMS 0
#define INCLUDE_GHOST_ELEMS 1
extern int MY_RANK,SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retval,iele,j,nel,k,ndof,p,nel,q,ndof)
#define MASK(j,k,p,q) VEC4(argd.profile,j,k,ndof,p,nel,q,ndof)
#define ICONE(j,k) VEC2(icone,j,k,nel)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "upload_vector"
#if 1
// New Fast PETSc matrix loading version (uses MatSetValues)

extern int any_A_LL_other_stop;

int Elemset::upload_vector(int nel,int ndof,Dofmap *dofmap,
		  int options,arg_data &argd,int myrank,
		  int el_start,int el_last,int iter_mode,
		  int klocc,int kdofc) {

  int iele,kloc,node,kdof,locdof,lloc,nodel,ldof,locdofl,ierr,
    load_vec,load_mat,load_mat_col,comp_prof,iele_here,
    pfmat;

  double *retval = argd.retval;

  int neq;
  //  sp_entry *spe, *spel;
  double val,vall;
  //  ierr = VecGetSize(vec,&neq); CHKERRQ(ierr); 
  neq = dofmap->neq;
  row_t::iterator entryc;
  
  // static Darray *row, *rowc;
  row_t row,rowc;
  IdMapRow row_v,rowc_v;
  IdMapEntry *entry_v,*entryc_v;

  load_vec = (options & (UPLOAD_VECTOR | UPLOAD_VECTOR_LOCST));
  load_mat = (options & (UPLOAD_MATRIX | UPLOAD_PROFILE));
  load_mat_col = (options & IS_FDJ);
  pfmat = (options & PFMAT);

  // In order to compute profiles (comp_prof==1) we have to do all the
  // work in all th processors
  comp_prof= (options & UPLOAD_PROFILE); // to be defined later

  InsertMode mode = ADD_VALUES;
  if (options & UPLOAD_VECTOR_LOCST) mode = INSERT_VALUES;
  
  dvector<int> row_map, row_map2, indxr, lnodr, dofr, indxc, lnodc, dofc;
  dvector<double> coefr, coefc, values;
  dvector<int> row_mask, col_mask;
  // These flag whether a given (lnod,dof) is
  // being loaded or not. By row or by column. 
  row_mask.mono(nel*ndof).reshape(2,nel,ndof).set(0);
  col_mask.mono(nel*ndof).reshape(2,nel,ndof).set(0);

  // This is an initial size (usually should be enough).
  // Anyway the vector grows if needed. 
  int rsize = 2*nel*ndof;
  // List of row eqs. to be loaded
  indxr.mono(rsize);
  // For each row eq. the local node index ...
  lnodr.mono(rsize);
  // row dof indx ...
  dofr.mono(rsize);
  // and coef (usually 1.)
  coefr.mono(rsize);

  int csize = rsize;
  // Same for columns
  indxc.mono(csize);
  // For each row eq. the local node index ...
  lnodc.mono(csize);
  // row dof indx ...
  dofc.mono(csize);
  // and coef (usually 1.)
  coefc.mono(csize);
  
  // This stores the values to be loaded in the matrix
  values.set_chunk_size((2*nel*ndof)*(2*nel*ndof));

  int nen = ndof*nel;
  int nen2 = nen*nen;

  int *indxrp = indxr.buff();
  int *lnodrp = lnodr.buff();
  int *dofrp = dofr.buff();
  double *coefrp = coefr.buff();

  int *indxcp = indxc.buff();
  int *lnodcp = lnodc.buff();
  int *dofcp = dofc.buff();
  double *coefcp = coefc.buff();

  if (load_mat) {
    // Compute row and column masks.
    for (kloc=0; kloc<nel; kloc++) {
      for (kdof=0; kdof<ndof; kdof++) {
	for (lloc=0; lloc<nel; lloc++) {
	  for (ldof=0; ldof<ndof; ldof++) {
	    if (MASK(kloc,kdof,lloc,ldof)) {
	      row_mask.e(kloc,kdof)=1;
	      col_mask.e(lloc,ldof)=1;
	    }
	  }
	}
      }
    }

#if 0
    for (kloc=0; kloc<nel; kloc++) {
      for (kdof=0; kdof<ndof; kdof++) {
	printf("(lnod=%d,dof=%d) row/col mask: %d %d\n",
	       kloc,kdof,row_mask.e(kloc,kdof),col_mask.e(kloc,kdof));
      }
    }
#endif
  }

  iele_here=-1;
  for (iele=el_start; iele<=el_last; iele++) {
    if (!compute_this_elem(iele,this,myrank,iter_mode)) continue;
    iele_here++;

    const int *dofv;
    const double *coefv;
    double coef;
    int n;
    if (load_mat) {
      // Build row map
      int jr=0,nr,jc=0,nc;
      int all_one_coef = 1;
      int *rm = row_mask.buff();
      int *cm = col_mask.buff();
      for (kloc=0; kloc<nel; kloc++) {
	node = ICONE(iele,kloc);
	for (kdof=0; kdof<ndof; kdof++) {
	  if (!*rm && !*cm) continue;
#if 0
          // Slow
	  dofmap->get_row(node,kdof+1,row_v);
	  for (int ientry=0; ientry<row_v.size(); ientry++) {
	    entry_v = &row_v[ientry];
	    locdof = entry_v->j;
	    coef = entry_v->coef;
	  }
#else
	  // Fast 
	  dofmap->get_row(node,kdof+1,n,&dofv,&coefv);
	  for (int j=0; j<n; j++) {
	    locdof = dofv[j];
	    coef = coefv[j];
#endif
	    if (locdof>neq) continue; // ony load free nodes
	    if (*rm) {
	      if (jr==rsize) {
		// Make arrays grow if needed
		rsize = 2*rsize;
		indxr.resize(rsize);
		lnodr.resize(rsize);
		dofr.resize(rsize);
		coefr.resize(rsize);

		indxr.defrag();
		lnodr.defrag();
		dofr.defrag();
		coefr.defrag();

		indxrp = indxr.buff();
		lnodrp = lnodr.buff();
		dofrp = dofr.buff();
		coefrp = coefr.buff();
	      }
	      indxrp[jr] = locdof-1;
	      lnodrp[jr] = kloc;
	      dofrp[jr] = kdof;
	      coefrp[jr] = coef;
	      if (coef!=1.0) all_one_coef = 0;
	      jr++;
	    }
	    if (*cm) {
	      if (jc==rsize) {
		// Make arrays grow if needed
		csize = 2*csize;
		indxc.resize(csize);
		lnodc.resize(csize);
		dofc.resize(csize);
		coefc.resize(csize);

		indxc.defrag();
		lnodc.defrag();
		dofc.defrag();
		coefc.defrag();
		indxcp = indxc.buff();
		lnodcp = lnodc.buff();
		dofcp = dofc.buff();
		coefcp = coefc.buff();
	      }
	      indxcp[jc] = locdof-1;
	      lnodcp[jc] = kloc;
	      dofcp[jc] = kdof;
	      coefcp[jc] = coef;
	      if (coef!=1.0) all_one_coef = 0;
	      jc++;
	    }
	  }
	  rm++;
	  cm++;
	}
      }
      nr = jr;
      nc = jc;

#if 0
      printf("iele %d\n",iele);
      for (jr=0; jr<nr; jr++) {
	printf("jr %d, lnod %d, dof %d, coef %f\n",jr,
	       lnodr.e(jr),dofr.e(jr),coefr.e(jr));
      }
      for (jc=0; jc<nc; jc++) {
	printf("jc %d, lnod %d, dof %d, coef %f\n",jc,
	       lnodc.e(jc),dofc.e(jc),coefc.e(jc));
      }
#endif

      if (comp_prof) {
	for (jr=0; jr<nr; jr++) {
	  locdof = indxr.e(jr);
	  for (jc=0; jc<nc; jc++) {
	    locdofl = indxc.e(jc);
	    // printf("set prof: %d, %d\n",locdof-1,locdofl-1);
	    if (pfmat) {
	      argd.pfA->set_profile(locdof,locdofl);
	    } else {
	      node_insert(argd.da,locdof,locdofl);
	      node_insert(argd.da,locdofl,locdof); 
	    }
	  }
	}
      } else {
	values.a_resize(2,nr,nc);
	values.defrag();
#if 1
	if (!all_one_coef) {
	  // New fast version (using C vector access)
	  double *w = values.buff();
	  for (jr=0; jr<nr; jr++) {
	    int lnodr1 = lnodrp[jr];
	    int dofr1 = dofrp[jr];
	    double coefr1 = coefrp[jr];
	    double *rtvm_row = &RETVALMAT(iele_here,lnodr1,dofr1,0,0);
	    for (jc=0; jc<nc; jc++) {
	      int lnodc1 = lnodcp[jc];
	      int dofc1 = dofcp[jc];
	      double coefc1 = coefcp[jc];
	      // if (!MASK(lnodr1,dofr1,lnodc1,dofc1)) continue;
	      // *w++ = coefr1*coefc1*RETVALMAT(iele_here,lnodr1,dofr1,lnodc1,dofc1);
	      *w++ = coefr1 * coefc1 * *(rtvm_row + lnodc1*ndof + dofc1);
	    }      
	  }
	} else {
	  // This block takes 0.9 secs in the
	  // cubiv cavity (cubcav) 30x30x30
	  // New fast version (using C vector access)
	  double *w = values.buff();
	  int *lnodr1 = lnodrp;
	  int *lnodr1_end = lnodrp+nr;
	  int *dofr1 = dofrp;
	  double *rtvm_iele = retval + iele_here*nen2;
	  while (lnodr1<lnodr1_end) {
	    double *rtvm_row = rtvm_iele 
	      + nen * ((*lnodr1++) * ndof + (*dofr1++));
	    int *lnodc1 = lnodcp;
	    int *lnodc1_end = lnodc1 + nc;
	    int *dofc1 = dofcp;
	    while (lnodc1 < lnodc1_end) {
	      // if (!MASK(lnodr1,dofr1,lnodc1,dofc1)) continue;
	      // *w++ = coefr1*coefc1*RETVALMAT(iele_here,lnodr1,dofr1,lnodc1,dofc1);
	      *w++ = *(rtvm_row + (*lnodc1++)*ndof + (*dofc1++));
	    }      
	  }
	}
#else
	// Old slow version (using dvector<int>.e())
	for (jr=0; jr<nr; jr++) {
	  int lnodr1 = lnodr.e(jr);
	  int dofr1 = dofr.e(jr);
	  double coefr1 = coefr.e(jr);
	  for (jc=0; jc<nc; jc++) {
	    int lnodc1 = lnodc.e(jc);
	    int dofc1 = dofc.e(jc);
	    double coefc1 = coefc.e(jc);
	    if (!MASK(lnodr1,dofr1,lnodc1,dofc1)) continue;
	    values.e(jr,jc) = coefr1*coefc1
	      *RETVALMAT(iele_here,lnodr1,dofr1,lnodc1,dofc1);
	  }      
	}
#endif
	if (pfmat) {
  	  ierr = argd.pfA->set_values(nr,indxr.buff(),nc,indxc.buff(),
  				      values.buff(),ADD_VALUES); 
  	  CHKERRQ(ierr); 
	} else {
	  ierr = MatSetValues(*argd.A,nr,indxr.buff(),nc,indxc.buff(),
			      values.buff(),ADD_VALUES);
	  CHKERRQ(ierr); 
	}
      }
    }
 
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(iele,kloc);
      for (kdof=0; kdof<ndof; kdof++) {
	dofmap->get_row(node,kdof+1,row_v);

	for (int ientry=0; ientry<row_v.size(); ientry++) {
	  entry_v = &row_v[ientry];
	  locdof = entry_v->j;
	  if (locdof>neq) continue; // ony load free nodes
	  val = (entry_v->coef) * RETVAL(iele_here,kloc,kdof);

	  // load local values on global vector
	  if (load_vec) {
	    VecSetValue(*(argd.x),locdof-1,val,mode);
	  }

	  // load finite difference jacobian computed by differences
	  if (load_mat_col) {
	    nodel = ICONE(iele,klocc);
	    vall = RETVAL(iele_here,kloc,kdof);
	    dofmap->get_row(nodel,kdofc+1,rowc_v);

	    for (int ientryc=0; ientryc<rowc_v.size(); ientryc++) {
	      entryc_v = &rowc_v[ientryc];
	      locdofl = entryc_v->j;
	      if (locdofl>neq) continue;
	      val = entry_v->coef * entryc_v->coef * vall;

	      int kd=locdof-1,kdl=locdofl-1;
	      if (MASK(kloc,kdof,klocc,kdofc)) {
		if (comp_prof) {
		  if (pfmat) {
		    argd.pfA->set_profile(kd,kdl);
		  } else {
		    node_insert(argd.da,kd,kdl);
		    node_insert(argd.da,kdl,kd); // be sure that
		    // profile
		  }
				// is symmetric
		} else {
		  if (pfmat) {
		    ierr = argd.pfA->set_value(kd,kdl,val,ADD_VALUES); 
		    CHKERRQ(ierr); 
		  } else {
		    MatSetValue(*argd.A,kd,kdl,val,ADD_VALUES); 
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
#if 0
  if (load_mat) {
    PetscSynchronizedFlush(PETSC_COMM_WORLD); 
    int any_A_LL_other_stop_all;
    MPI_Allreduce(&any_A_LL_other_stop,
		  &any_A_LL_other_stop_all,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
    if (any_A_LL_other_stop_all) {
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			      "[%d] any_A_LL_other: %d\n",MY_RANK,any_A_LL_other_stop);
      PetscSynchronizedFlush(PETSC_COMM_WORLD); 
      PetscFinalize();
      exit(0);
    }
  }
#endif
}
#else
// Old slow PETSc matrix loading version (uses MatSetValue)
int Elemset::upload_vector(int nel,int ndof,Dofmap *dofmap,
		  int options,arg_data &argd,int myrank,
		  int el_start,int el_last,int iter_mode,
		  int klocc,int kdofc) {

  //  printf("Elemset::upload_vector: old single element version\n");
  int iele,kloc,node,kdof,locdof,lloc,nodel,ldof,locdofl,ierr,
    load_vec,load_mat,load_mat_col,comp_prof,iele_here,
    pfmat;

  double *retval = argd.retval;

  int neq;
  //  sp_entry *spe, *spel;
  double val,vall;
  //  ierr = VecGetSize(vec,&neq); CHKERRQ(ierr); 
  neq = dofmap->neq;
  row_t::iterator entryc;
  
  // static Darray *row, *rowc;
  row_t row,rowc;
  IdMapRow row_v,rowc_v;
  IdMapEntry *entry_v,*entryc_v;

  load_vec = (options & (UPLOAD_VECTOR | UPLOAD_VECTOR_LOCST));
  load_mat = (options & (UPLOAD_MATRIX | UPLOAD_PROFILE));
  load_mat_col = (options & IS_FDJ);
  pfmat = (options & PFMAT);

  // In order to compute profiles (comp_prof==1) we have to do all the
  // work in all th processors
  comp_prof= (options & UPLOAD_PROFILE); // to be defined later

  InsertMode mode = ADD_VALUES;
  if (options & UPLOAD_VECTOR_LOCST) mode = INSERT_VALUES;

  iele_here=-1;
  for (iele=el_start; iele<=el_last; iele++) {
    if (!compute_this_elem(iele,this,myrank,iter_mode)) continue;
    iele_here++;
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(iele,kloc);
      for (kdof=0; kdof<ndof; kdof++) {
	dofmap->get_row(node,kdof+1,row_v);

	for (int ientry=0; ientry<row_v.size(); ientry++) {
	  entry_v = &row_v[ientry];
	  locdof = entry_v->j;
	  if (locdof>neq) continue; // ony load free nodes
	  val = (entry_v->coef) * RETVAL(iele_here,kloc,kdof);

	  // load local values on global vector
	  if (load_vec) {
	    VecSetValue(*(argd.x),locdof-1,val,mode);
	  }

	  // load finite difference jacobian computed by differences
	  if (load_mat_col) {
	    nodel = ICONE(iele,klocc);
	    vall = RETVAL(iele_here,kloc,kdof);
	    dofmap->get_row(nodel,kdofc+1,rowc_v);

	    for (int ientryc=0; ientryc<rowc_v.size(); ientryc++) {
	      entryc_v = &rowc_v[ientryc];
	      locdofl = entryc_v->j;
	      if (locdofl>neq) continue;
	      val = entry_v->coef * entryc_v->coef * vall;

	      int kd=locdof-1,kdl=locdofl-1;
	      if (MASK(kloc,kdof,klocc,kdofc)) {
		if (comp_prof) {
		  if (pfmat) {
		    argd.pfA->set_profile(kd,kdl);
		  } else {
		    node_insert(argd.da,kd,kdl);
		    node_insert(argd.da,kdl,kd); // be sure that
		    // profile
		  }
				// is symmetric
		} else {
		  if (pfmat) {
		    ierr = argd.pfA->set_value(kd,kdl,val,ADD_VALUES); 
		    CHKERRQ(ierr); 
		  } else {
		    MatSetValue(*argd.A,kd,kdl,val,ADD_VALUES); 
		  }
		}
	      }
	    }
	  }

	  // load local values on global matrix
	  if (!load_mat) continue;
	  for (lloc=0; lloc<nel; lloc++) {
	    nodel = ICONE(iele,lloc);
	    for (ldof=0; ldof<ndof; ldof++) {
	      vall = RETVALMAT(iele_here,kloc,kdof,lloc,ldof);
	      dofmap->get_row(nodel,ldof+1,rowc_v);

	      for (int ientryc=0; ientryc<rowc_v.size(); ientryc++) {
		entryc_v = &rowc_v[ientryc];
		locdofl = entryc_v->j;
		if (locdofl>neq) continue; // only load for free dof's
		
		val = (entry_v->coef) * (entryc_v->coef) * vall;
		if (MASK(kloc,kdof,lloc,ldof)) {
		  if (comp_prof) {
		    // be sure that profile is symmetric
		    if (pfmat) {
		      argd.pfA->set_profile(locdof-1,locdofl-1);
		    } else {
		      node_insert(argd.da,locdof-1,locdofl-1);
		      node_insert(argd.da,locdofl-1,locdof-1); 
		    }
		  } else {
		    // printf("(%d,%d) -> %f\n",locdof,locdofl,val);
		    if (pfmat) {
		      ierr = argd.pfA->set_value(locdof-1,locdofl-1,val,ADD_VALUES); 
		      CHKERRQ(ierr); 
		    } else {
		      MatSetValue(*argd.A,locdof-1,locdofl-1,val,ADD_VALUES);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}
#endif
