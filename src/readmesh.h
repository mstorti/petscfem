// -*- mode: c++ -*-

/*
  This file belongs to the PETSc - FEM package, a library and
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
 
#ifndef READ_MESH_H
#define READ_MESH_H

/** @name read\_mesh package */
//@{
/** This struct contains the info for the `props' hash table which
    gives the position in the per-element props table for a given
    text. 
    @author M. Storti
    @param position position (column) in the table
    @param width number of scalars
*/ 
struct props_hash_entry {
public:
  int position,width;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Bless elemset with the given type (derived class). 

    @author M. Storti
    @param type type of elemset
    @param elemset to be blessed
*/ 
void bless_elemset(char *type,Elemset *& elemset);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** reads mesh and returns the corresponding mesh and dofmap. 
    @author M. Storti

    @param mesh (output) the mesh that has been read
    @param fcase (input) name of the file to be read
    @param dofmap (output) the dofmap read
    @param neq (output) total number of unknowns. 
    @param size (input) number of processors running
    @param myrank (input) this processor number (base 0)
    @param x (output) reduced (MPI) vector prototype
    @param xseq (output) reduced (sequential)  vector prototype
*/ 
int read_mesh(Mesh *& mesh,char *fcase,Dofmap *& dofmap,
	      int & neq,int size,int myrank);

//@}

/** @name read\_mesh package */

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** For a given state vector, prints the state for some nodes.
    @author M. Storti
    @param (input) filename file where to write the vector. May contain
    relative directories. 
    @param x (input) PETSc MPI vector to be written 
    @param dofmap (input) corresponding dofmap 
    @param time_data (input, def=NULL) an external parameter in order to compute
    external boundary conditions, etc...
    @param node_list (input) set of nodes to print. 
    @param append (input, def=0) appending mode (append if
    `append==0') 
*/ 
int print_some(const char *filename,Vec x,Dofmap *dofmap,
	       set<int> node_list,const TimeData *time_data=NULL);

/** Prints a vector to a file
    @author M. Storti
    @param filenamepat (input) pattern to generate the filename (must contain a \%d)
    @param x (input) the vector to be printed
    @param dofmap (input) the Dofmap of the problem
    @param time_data (input) an external parameter in order to compute
    external boundary conditions, etc...
    @param j (input) the time step
    @param nsave (input) save each nsave time steps
    @param nrec (input) save nrec records per file
    @param nfile (input) rotate through `nfile' files
*/
void print_vector_rota(const char *filenamepat,const Vec x,const
		       Dofmap *dofmap,const TimeData *time_data,
		       const int j,const int nsave,const int nrec,
		       const int nfile);

#endif
 
