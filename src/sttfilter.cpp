/* $Id: sttfilter.cpp,v 1.2 2001/01/19 12:46:32 mstorti Exp $ */

/*
  This file belongs to he PETSc - FEM package a library and
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

#include "sttfilter.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
#undef __FUNC__
#define __FUNC__ "Filter::Filter(Vec &,Mesh &)"
Filter::Filter(Vec &x,Mesh &mesh) {
  // Read values from Texthash
  //o Flag defining whether to define a average filter
  TGETOPTDEF_ND(mesh.global_options,int,"filter_av",0);

  if (filter_av) {
    ierr = VecDuplicate(x,xav); CHKERRA(ierr);
  }

}
#endif

void LowPass::update() {
  if (input->step() < step()) input->update(); 
  i_state.scale(1-alpha).axpy(alpha,input->state());
  Filter::update();
}

const State & Inlet::state() { 
  return *state_;
}

const State & LowPass::state() {
  return i_state;
}

const State & Mixer::state() {
  return i_state;
}

Mixer & Mixer::add_input(Filter &input,double g) {
  filter_l.push_back(&input);
  gain_l.push_back(g);
  return *this;
}

#if 0
void Mixer::update() {
  input->update(); 
  if (cmixer) cmixer->update();
  if (!cmixer) *i_state = 0;
  *i_state += g * input->state(); 
}
#endif

void Mixer::update() {
  i_state.set_cnst(0.);
  for (int j=0; j<filter_l.size(); j++) {
    Filter *f = filter_l[j];
    if (f->step() < step()) f->update(); 
    f->update();
    i_state.axpy(gain_l[j],f->state());
  }
  Filter::update();
  // time_step++;
}

