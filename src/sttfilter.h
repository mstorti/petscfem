// -*-mode: c++ -*-
/* $Id: sttfilter.h,v 1.1 2001/01/18 11:32:56 mstorti Exp $ */

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
 
#ifndef STTFILTER_H
#define STTFILTER_H

#include <cstdio>
#include <cmath>
#include <vector>

typedef double FilterState;

class Filter {
  int time_step;
 public:
  int step() {return time_step;}
  Filter(int ts=0) {time_step=ts;}
  virtual void update() {time_step++;};
  virtual const FilterState & state()=0;
};

class Inlet : public Filter {
  const FilterState *state_;
public:
  Inlet(FilterState &st) { state_ = &st;}
  const FilterState & state();
};

class LowPass : public Filter {
  // Input to the filter
  Filter *input;
  // Internal state
  FilterState i_state; 
  // Relaxation factor;
  double alpha;
public:
  LowPass(double alpha_,Filter &input_,FilterState &state_) : alpha(alpha_), input(&input_), 
    i_state(state_) {};
  void update();
  const FilterState & state();
};

class Mixer : public Filter {
  // Input to the filter
  vector<Filter *> filter_l;
  vector<double> gain_l;
  // Internal state
  FilterState i_state; 
public:
  Mixer(FilterState &st) {i_state = st;}
  Mixer & add_input(Filter &input,double g);
  void update();
  const FilterState & state();
  ~Mixer() {};
};

#endif

