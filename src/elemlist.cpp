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
  return ElementIterator(this,first,0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
ElementIterator ElementList::end(void) const {
  return ElementIterator(this,last,0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
ElementIterator & ElementIterator::operator++(void) {
  
  while (++rank_in_elemset < elemlist->last) {
    if (elemlist->elemset->epart[rank_in_elemset] == MY_RANK+1) break;
    if (elemlist->mode == INCLUDE_GHOST) {
      assert(1);
#if 0   // to be converted
      int is_ghost_elem = da_bsearch(elemset->ghost_elems,&iele,int_cmp,NULL);
      return (is_ghost_elem>=0);
#endif
    }
  }
  rank_in_chunk++;
  return *this;
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

