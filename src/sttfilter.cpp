/* $Id: sttfilter.cpp,v 1.3 2001/01/22 00:51:47 mstorti Exp $ */

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

void LowPass::update() {
  if (input->step() <= step()) input->update(); 
  i_state.scale(gamma).axpy(1-gamma,input->state());
  Filter::update();
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
    if (f->step() <= step()) f->update(); 
    f->update();
    i_state.axpy(gain_l[j],f->state());
  }
  Filter::update();
}
