//__INSERT_LICENSE__
//$Id: elmsupl.cpp,v 1.4 2003/08/28 20:42:39 mstorti Exp $

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
  // int rsize = 2*nel*ndof;
  int rsize = 1;
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

  // Compute row and column masks.
  for (kloc=0; kloc<nel; kloc++) {
    for (kdof=0; kdof<ndof; kdof++) {
      for (lloc=0; lloc<nel; lloc++) {
	for (ldof=0; ldof<ndof; ldof++) {
	  if (MASK(kloc,kdof,klocc,kdofc)) {
	    row_mask.e(kloc,kdof)=1;
	    col_mask.e(lloc,ldof)=1;
	  }
	}
      }
    }
  }

  for (kloc=0; kloc<nel; kloc++) {
    for (kdof=0; kdof<ndof; kdof++) {
      printf("(lnod=%d,dof=%d) row/col mask: %d %d\n",
	     kloc,kdof,row_mask.e(kloc,kdof),col_mask.e(kloc,kdof));
    }
  }

  iele_here=-1;
  for (iele=el_start; iele<=el_last; iele++) {
    if (!compute_this_elem(iele,this,myrank,iter_mode)) continue;
    iele_here++;

    // Build row map
    int jr=0,nr;
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(iele,kloc);
      for (kdof=0; kdof<ndof; kdof++) {
	if (!row_mask.e(kloc,kdof)) continue;
	dofmap->get_row(node,kdof+1,row_v);
	for (int ientry=0; ientry<row_v.size(); ientry++) {
	  entry_v = &row_v[ientry];
	  locdof = entry_v->j;
	  if (locdof>neq) continue; // ony load free nodes
	  if (jr==rsize) {
	    // Make arrays grow if needed
	    printf("resize %d -> %d\n",rsize,2*rsize);
	    rsize = 2*rsize;
	    indxr.resize(rsize);
	    lnodr.resize(rsize);
	    dofr.resize(rsize);
	    coefr.resize(rsize);
	  }
	  indxr.e(jr) = locdof;
	  lnodr.e(jr) = kloc;
	  dofr.e(jr) = kdof;
	  coefr.e(jr) = entry_v->coef;
	  jr++;
	}
      }
    }
    nr = jr;

#if 1
    printf("iele %d\n",iele);
    for (jr=0; jr<nr; jr++) {
      printf("jr %d, lnod %d, dof %d, coef %f\n",jr,
	     lnodr.e(jr),dofr.e(jr),coefr.e(jr));
    }
#endif
    PetscFinalize();
    exit(0);
 
    for (kloc=0; kloc<nel; kloc++) {
      for (kdof=0; kdof<ndof; kdof++) {
	dofmap->get_row(node,kdof+1,row_v);

	for (int ientry=0; ientry<row_v.size(); ientry++) {
	  entry_v = &row_v[ientry];
	  locdof = entry_v->j;
	  if (locdof>neq) continue; // ony load free nodes
	  val = (entry_v->coef) * RETVAL(iele_here,kloc,kdof);
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
		      ierr = 0;
		      // ierr = argd.pfA->set_value(locdof-1,locdofl-1,val,ADD_VALUES); 
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
#else
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
