// -*-mode: c++ -*-
/* $Id: sttfilter.h,v 1.3 2001/01/22 00:51:47 mstorti Exp $ */

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
#include "readmesh.h"

class State {
  Vec *vec;
  Time time;
public:
  State(Vec &vec,Time t) : vec(&vec), time(t) {};
  State(const State &v);
  ~State();
  State & set_time(const Time t) {time = t;}
  State & axpy(double gamma,const State &v);
  State & scale(double gamma);
  State & set_cnst(double a);
  State & inc(double dt) {time.inc(dt);}
  const Vec & v() const {return *vec;}
  operator Vec &() {return *vec;}
  const Time & t() const {return time;}
  const State & print_some(const char *filename,Dofmap *dofmap,
		     set<int> & node_list) const;
};

class Filter {
  int time_step;
 public:
  int step() {return time_step;}
  Filter(int ts=0) {time_step=ts;}
  virtual void update() {time_step++;};
  virtual const State & state() const = 0;
  virtual operator const State &() const {return state();}
};

class Inlet : public Filter {
  const State *state_;
public:
  Inlet(State &st) { state_ = &st;}
  const State & state() const {return *state_;}
  ~Inlet() {};
  // operator const State &() const {return *state_;}
};

class LowPass : public Filter {
  // Input to the filter
  Filter *input;
  // Internal state
  State i_state; 
  // Relaxation factor;
  double gamma;
public:
  LowPass(double gamma_,Filter &input_,State &state_) 
    : gamma(gamma_), input(&input_), i_state(state_) {};
  ~LowPass();
  void update();
  const State & state() const {return i_state;}
  // operator const State &() const {return i_state;};
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
  //operator const State &() const {return i_state;};
  const State & state();
  ~Mixer() {};
};

#endif
