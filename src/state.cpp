/*
  This file belongs to the PETSc - FEM package, a library and
  application suite oriented to the Finite Element Method based on PETSc. 
  Copyright (C) 1999, 2000, 2001  Mario Alberto Storti
  
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

#include "sttfilter.h"
 
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "State::axpy(double,const State &)"
State & State::axpy(double alpha,const State &v) {
  int ierr = VecAXPY(&alpha,*(v.vec),*(this->vec));
  return *this;
}

State::State(const State &v) : time(v.time) {
  int ierr = VecDuplicate(*vec,v.vec);
}
  
State::~State() {
  int ierr = VecDestroy(*vec); 
}
  
State & State::scale(double alpha) {
  double alpha_ = alpha;
  int ierr = VecScale(&alpha_,*vec);
  return *this;
}

State & State::set_cnst(double a) {
  double a_ = a;
  int ierr = VecSet(&a,*vec);
  return *this;
}
