// -*- mode: C++ -*-

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

#ifndef ARGLISTN_H
#define ARGLISTN_H

#include <vector>
#include <string>
#include <mat.h>
#include <vec.h>
#include "libretto.h"

class ArgEntry {
public:
  virtual void chunk_pre_operations(int chunk_size,int ndoft) {};
};

class OutArg : public ArgEntry {
public:
  virtual void chunk_pre_operations(int chunk_size,int ndoft);
private:
  double *retval;
}

class OutVectorArg : public ArgEntry {
public:
  OutVectorArg(Vec &x);
  ~OutVectorArg();
  void chunk_pre_operations(int chunk_size,int ndoft);
private:
  double *retval;
}

typedef vector<ArgEntry *> ArgList;

#endif
