/*__INSERT_LICENSE__*/
//$Id: elemlist.cpp,v 1.10 2001/05/30 03:58:50 mstorti Exp $

#include "fem.h"
#include <set>
#include "utils.h"
#include "getprop.h"
#include "elemset.h"
#include "idmap.h"
#include "dofmap.h"
#include "arglist.h"

// iteration modes
#define NOT_INCLUDE_GHOST_ELEMS 0
#define INCLUDE_GHOST_ELEMS 1
extern int MY_RANK,SIZE;

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
ElementList::set(const ElemsetIteratorMode mode_,const int
		 first_,const int last_) {
  mode = mode_;
  first=first_;
  last=last_;
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
ElementIterator ElementList::begin(void) const {
  ElementIterator tmp=ElementIterator(this,first,0);
  tmp.advance_to_valid();    
  return tmp;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
ElementIterator ElementList::end(void) const {
  return ElementIterator(this,last,0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
ElementIterator & ElementIterator::operator++(void) {
  rank_in_elemset++;
  advance_to_valid();
  rank_in_chunk++;
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int ElementIterator::is_valid(void) const {
#if 0
  if (elemlist->elemset->epart[rank_in_elemset] == MY_RANK+1) return 1;
  if (elemlist->mode == INCLUDE_GHOST) {
    int is_ghost_elem = da_bsearch(elemlist->elemset->ghost_elems,
				   &rank_in_elemset,int_cmp,NULL);
    return (is_ghost_elem>=0);
  }
  return 0;
#else
  return compute_this_elem(rank_in_elemset,elemlist->elemset,
			   MY_RANK,elemlist->mode==INCLUDE_GHOST);
#endif
  
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void ElementIterator::advance_to_valid(void) {
  while (rank_in_elemset < elemlist->last && !is_valid()) {
    ++rank_in_elemset;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
ElementIterator ElementIterator::operator++(int) {
  ElementIterator tmp=*this;
  ++(*this);
  return tmp;
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int ElementList::iterator::not_at_end(void) const {
  return global<last;
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int ElementIterator::operator!=(const ElementIterator &other) const {
  return rank_in_elemset!=other.rank_in_elemset;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int 
NewElemset::assemble(arg_data_list &arg_datav,Nodedata *nodedata,Dofmap *dofmap,
		     const char *jobinfo,int myrank,
		     int el_start,int el_last,int iter_mode,
		     const TimeData *time_data) {

  ElementList elemlist(this,el_start,el_last+1,
		       (iter_mode==INCLUDE_GHOST_ELEMS ? 
			INCLUDE_GHOST : DO_NOT_INCLUDE_GHOST));
  new_assemble(arg_datav,nodedata,dofmap,
	       jobinfo,elemlist,time_data);
  return 0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "element_connect"
void Elemset::element_connect(const ElementIterator &element,
			      int *connect) const {
  int k,ielh;
  element.position(k,ielh);
  for (int kloc=0; kloc<nel; kloc++) {
    connect[kloc] = icone[nel*k+kloc];
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Elemset::element_props"
double *
Elemset::element_props(const ElementIterator &element) const {
  int k,ielh;
  element.position(k,ielh);
  return elemprops+k*nelprops;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "ElementIterator::props"
double *
ElementIterator::props() {
  return elemlist->elemset->element_props(*this);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "element_node_data"
void 
Elemset::element_node_data(const ElementIterator &element,
			   const Nodedata *nodedata,
			   double *xloc,double *Hloc) const {
#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nodedata->nu)
#define XLOC(j,k) VEC2(xloc,j,k,nodedata->ndim)
#define HLOC(j,k) VEC2(Hloc,j,k,nH)
  element_connect(element,elem_conne);
  int nH = nodedata->nu - nodedata->ndim;
  for (int kloc=0; kloc<nel; kloc++) {
    int node = elem_conne[kloc];
    for (int jdim=0; jdim<nodedata->ndim; jdim++)
      XLOC(kloc,jdim) = NODEDATA(node-1,jdim);
    for (int ih=0; ih<nH; ih++) 
      HLOC(kloc,ih) = NODEDATA(node-1,nodedata->ndim+ih);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "ElementIterator::node_data"
void ElementIterator::node_data(const Nodedata *nodedata,
				double *xloc,double *Hloc) {
  elemlist->elemset->element_node_data(*this,nodedata,xloc,Hloc);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Elemset::element_vector_values"
const double *
Elemset::element_vector_values(const ElementIterator &element,
			       arg_data &ad) const {
  int k,ielh;
  element.position(k,ielh);
  return (ad.locst)+ielh*nel*ndof;
}
  
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "ElementIterator::vector_values"
const double *
ElementIterator::vector_values(arg_data &ad) const {
  return elemlist->elemset->element_vector_values(*this,ad);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Elemset::element_ret_vector_values"
double *
Elemset::element_ret_vector_values(const ElementIterator &element,
				   arg_data &ad) const {
  int k,ielh;
  element.position(k,ielh);
  return (ad.retval)+ielh*nel*ndof;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "ElementIterator::ret_vector_values"
double *
ElementIterator::ret_vector_values(arg_data &ad) const {
  return elemlist->elemset->element_ret_vector_values(*this,ad);
}
  
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Elemset::element_ret_fdj_values"
double *
Elemset::element_ret_fdj_values(const ElementIterator &element,
				arg_data &ad) const {
  int k,ielh;
  element.position(k,ielh);
  return (ad.retval)+ielh*nel*ndof;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "ret_fdj_values"
double *
ElementIterator::ret_fdj_values(arg_data &ad) const {
  return elemlist->elemset->element_ret_fdj_values(*this,ad);
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Elemset::element_ret_mat_values"
double *
Elemset::element_ret_mat_values(const ElementIterator &element,
				arg_data &ad) const {
  int k,ielh;
  element.position(k,ielh);
  int tmp=nel*ndof;
  return (ad.retval)+ielh*tmp*tmp;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "ret_mat_values"
double *
ElementIterator::ret_mat_values(arg_data &ad) const {
  return elemlist->elemset->element_ret_mat_values(*this,ad);
}  

