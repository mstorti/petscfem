// -*- mode: c++ -*-
 
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
 
#ifndef LAPLA_H
#define LAPLA_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class lapla : public Elemset {
public:
  int assemble(arg_data_list &arg_data_v,Nodedata *nodedata,Dofmap *dofmap,
	       char *jobinfo,int myrank,
	       int el_start,int el_last,int iter_mode,
	       const TimeData *time_data);
};
#endif

