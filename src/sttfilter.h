// -*-mode: c++ -*-
/* $Id: sttfilter.h,v 1.5 2001/01/23 13:26:35 mstorti Exp $ */

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
#include "utils.h"

/** States are a state vectr plus a time, so that they can be passed
    to a dofmap and we can have all the nodal values from it. In a
    future we can have a generic StateFilter (double arrays) class so
    that Filters can get their values from them. 
    @author M. Storti
*/
class State {
  /// The state vector
  Vec *vec;
  /// The time corresponding to this state vector. Values on Dirichlet
  /// boundaries with time dependednt are obtained using these 
  Time time;
public:
  /// Constructor from state vector and time
  State(Vec &vec,Time t) : vec(&vec), time(t) {};
  /// Constructor from another state
  State(const State &v);
  /// Destructor
  ~State();

  /**@name Operations on the time part */
  //@{
  /// Changes the time of the state
  State & set_time(const Time t) {time = t;}
  /// Increments the time part
  State & inc(double dt) {time.inc(dt);}
  /// Const access to the time part
  const Time & t() const {return time;}
  //@}

  /**@name Operations on the state vector part */
  //@{ 
  /// axpy operation on the state vector part: \verb+ *this += gamma * v+
  State & axpy(double gamma,const State &v);
  /// scales state vector
  State & scale(double gamma);
  /// sets to a constant
  State & set_cnst(double a);
  /// Const acces to the vector part
  const Vec & v() const {return *vec;}
  /// Converts to a Petsc vector
  operator Vec &() {return *vec;}
  //@}

  /**@name Utilities */
  //@{ 
  /// Print some node/field combinations
  const State & print_some(const char *filename,Dofmap *dofmap,
		     set<int> & node_list) const;
  /// Print the whole state vector
  const State & print() const;
  /// Prints only the first n items of the local part of the state
  /// vector, mainly for debugging. 
  const State & print(int n) const;
};

/** This is the basic, virtual filter class.
    @author M. Storti
*/
class Filter {
  /// The actual time step
  int time_step;
 public:
  /// Returns the time_step
  int step() {return time_step;}
  /// Default constructor
  Filter(int ts=0) {time_step=ts;}
  /// Updates the filter from its inputs
  virtual void update(Time time) {time_step++;};
  /// Returns the state of the filter
  virtual const State & state() const = 0;
  /// Converts to state
  virtual operator const State &() const {return state();}
};

/** This is the basic filter that connects a external state to the
    filter chain. It has no internal state. 
    @author M. Storti
*/
class Inlet : public Filter {
  /// The external state from where the signal is received
  const State *state_;
public:
  /// Constructor from the external signal
  Inlet(State &st) { state_ = &st;}
  /// Instantiation of the state function
  const State & state() const {return *state_;}
  /// Destructor
  ~Inlet() {};
};

/** This is a recursive filter of the form 
    \TEX{$\hat u_{j+1} = \gamma \hat u_j + (1-\gamma) u_{j+1}$}.
    One usually enters \TEX{$\alpha$} such that \TEX{$\gamma = \exp{-\alpha\Delta t}$}. 
    @author M. Storti
*/
class LowPass : public Filter {
  // Input to the filter
  Filter *input;
  // Internal state
  State i_state; 
  // Relaxation factor;
  double gamma;
public:
  /// Constructor 
  LowPass(double gamma_,Filter &input_,State &state_) 
    : gamma(gamma_), input(&input_), i_state(state_) {};
  /// Destructor
  ~LowPass() {};
  /// Updates the filter (and its input)
  void update(Time time);
  /// Allows access to the internal state
  State & state() {return i_state;}
  /// Allows const access to the internal state
  const State & state() const {return i_state;}
};

/** The Mixer class should not have an internal state and is a 
    linear combination of its inputs. 
    @author M. Storti
*/
class Mixer : public Filter {
  /// List of inputs
  vector<Filter *> filter_l;
  /** List of gain factors. The output 
      is #output = sum_j {gain_l[j] * filter_l[j].state}#
  */
  vector<double> gain_l;
  /// Internal state. 
  State i_state; 
public:
  /// Constructor from a typical state
  Mixer(State &st) : i_state(st) {};
  /// Destructor
  ~Mixer() {};
  /// Adds an input to the list of inputs
  Mixer & add_input(Filter &input,double g);
  /// Updates the state and all its inputs. 
  void update(Time time);
  /// Gives access to the internal state
  State & state() {return i_state;}
  /// Gives const access to the internal state
  const State & state() const {return i_state;}
};

/** Family of LowPass filters.
    Contains typically several chains of low pass 
    filters at different orders. 
 */
class LPFilterGroup {
  /// Input to the filter chain
  Inlet input;
  /// List of lists of filters
  vector <vector < LowPass *> > lp_filters;
  /// List of relaxation factors
  vector<double> gamma_v;
  /// Orders (length) of each chain
  vector<int> n_v;
  /// Number of different chains
  int nalpha;
public:
  /// Constructor from the pointer to the general options
  LPFilterGroup(TextHashTable *thash,State &x,double Dt);
  /// Destructor
  ~LPFilterGroup();
  /// Updates of the filters in the chain
  void update(Time time);
  /// Returns the k-th filter in the j-th chain
  Filter & filter(int j, int k=1);
  /// Prints filtered states
  void print_some(const char *filename,Dofmap *dofmap,
		  set<int> & node_list);
};

#endif
