//__INSERT_LICENSE__
//$Id: elmsupl.cpp,v 1.1 2003/08/28 03:05:46 mstorti Exp $

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
  
  dvector<int> row_map, row_map2;
  dvector<double> coef2;
  row_map.mono(nel*ndof);
  row_map2.mono(nel*ndof);

  iele_here=-1;
  for (iele=el_start; iele<=el_last; iele++) {
    if (!compute_this_elem(iele,this,myrank,iter_mode)) continue;
    iele_here++;
    assert(!load_mat_col);

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
