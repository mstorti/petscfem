/*
  This file belongs to the PETSc - FEM package, a library and
  application suite oriented to the Finite Element Method based on PETSc. 
  Copyright (C) 1999, 2000  Mario Alberto Storti
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.

*/

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
		     char *jobinfo,int myrank,
		     int el_start,int el_last,int iter_mode,
		     const TimeData *time_data) {

  ElementList elemlist(this,el_start,el_last+1,
		       (iter_mode==INCLUDE_GHOST_ELEMS ? 
			INCLUDE_GHOST : DO_NOT_INCLUDE_GHOST));
  new_assemble(arg_datav,nodedata,dofmap,
	       jobinfo,elemlist,time_data);
  return 0;
};

