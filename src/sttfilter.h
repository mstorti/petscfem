// -*-mode: c++ -*-
/* $Id: sttfilter.h,v 1.2 2001/01/19 12:46:32 mstorti Exp $ */

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

#include <vec.h>

#include "fem.h"
#include "dofmap.h"

#if 0

typedef double State;

#else

class State {
  Vec *vec;
  Time time;
public:
  State(Vec &vec,Time t) : vec(&vec), time(t) {};
  State(const State &v);
  ~State();
  State & set_time(const Time t) {time = t;}
  State & axpy(double alpha,const State &v);
  State & scale(double alpha);
  State & set_cnst(double a);
  State & inc(double dt) {time.inc(dt);}
};

#endif

class Filter {
  int time_step;
 public:
  int step() {return time_step;}
  Filter(int ts=0) {time_step=ts;}
  virtual void update() {time_step++;};
  virtual const State & state()=0;
};

class Inlet : public Filter {
  const State *state_;
public:
  Inlet(State &st) { state_ = &st;}
  const State & state();
};

class LowPass : public Filter {
  // Input to the filter
  Filter *input;
  // Internal state
  State i_state; 
  // Relaxation factor;
  double alpha;
public:
  LowPass(double alpha_,Filter &input_,State &state_) : alpha(alpha_), input(&input_), 
    i_state(state_) {};
  void update();
  const State & state();
};

class Mixer : public Filter {
  // Input to the filter
  vector<Filter *> filter_l;
  vector<double> gain_l;
  // Internal state
  State i_state; 
public:
  Mixer(State &st) : i_state(st) {};
  Mixer & add_input(Filter &input,double g);
  void update();
  const State & state();
  ~Mixer() {};
};

#endif

