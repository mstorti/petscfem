/* $Id: sttfilter.cpp,v 1.1 2001/01/18 11:32:56 mstorti Exp $ */

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
  i_state = (1-alpha) * i_state + alpha * input->state();
  Filter::update();
  //  time_step++;
}

const FilterState & Inlet::state() { 
  return *state_;
}

const FilterState & LowPass::state() {
  return i_state;
}

const FilterState & Mixer::state() {
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
  i_state=0;
  for (int j=0; j<filter_l.size(); j++) {
    Filter *f = filter_l[j];
    if (f->step() < step()) f->update(); 
    f->update();
    i_state += gain_l[j] * f->state();
  }
  Filter::update();
  // time_step++;
}

int main() {
  FilterState x=0;
  Inlet i(x);
  LowPass f(.1,i,x),ff(.1,f,x);
  Mixer m(x);
  m.add_input(f,1.).add_input(ff,2.);
  
  double Dt=0.1,t=0,omega=1.;

  FILE *fid = fopen("filter.out","w");
  for (int j=0; j<100; j++) {
    t += Dt;
    x = sin(omega*t);
    fprintf(fid,"%f %f %f \n",f.state(),ff.state(),m.state());
    m.update();
  }
  fclose(fid);
}
