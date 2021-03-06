// -*- mode: c++ -*-
/*__INSERT_LICENSE__*/
//$Id: readmesh.h,v 1.12 2004/01/26 20:22:34 mstorti Exp $
 
#ifndef READ_MESH_H
#define READ_MESH_H

#include "sttfilter.h"

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
/** Bless elemset with the given type (derived class) for generic
    classes (not belonging to a specific application). 
    @author M. Storti
    @param type type of elemset
    @param elemset to be blessed
*/ 
void bless_elemset0(char *type,Elemset *& elemset);

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
    @param node_list (input) set of nodes to print. 
    @param time_data (input, def=NULL) an external parameter in order to compute
    external boundary conditions, etc...
*/ 
int print_some(const char *filename,const Vec x,Dofmap *dofmap,
	       set<int> node_list,const TimeData *time_data=NULL);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** For a given state vector, prints the state for some nodes.
    @author M. Storti
    @param (input) filename file where to write the vector. May contain
    relative directories. 
    @param s (input) State vector
    @param dofmap (input) corresponding dofmap 
    @param node_list (input) set of nodes to print. 
*/ 
#if 0
int print_some(const char *filename,const State &s,Dofmap *dofmap,
	       set<int> & node_list);
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Initializes the print\_some saving mechanism. 
    Basically reads the nodes to be printed at each time step. 
    @author M. Storti
    @param thash (input) the text hash from where to get properties. 
    @param print_some_file (input) the file where to read the list of
    nodes
    @param save_file_some (input) the file where to write the values
    at those nodes (this is not used currently)
    @param node_list (output) the list of nodes 
    @param save_file_some_append (input) Access mode to the ``some''
    file. If 0 rewind file. If 1 append to previous results.  
*/
int print_some_file_init(TextHashTable *thash,
			 const char *print_some_file,
			 const char *save_file_some,set<int> &node_list,
			 int save_file_some_append);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Reads a double value token from a string. 
    @author M. Storti
    @param rflag (input/output) First time pass initially 0 and then
    sets to 1.  
    @param line (input) the string to be read
    @param val (input) the double value to be read
    @return boolean value indicating wether a value has been read or
    dnot. 
*/ 
int readval(int &rflag,char *line,double &val);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Reads an int value token from a string. 
    @author M. Storti
    @param rflag (input/output) First time pass initially 0 and then
    sets to 1.  
    @param line (input) the string to be read
    @param val (output) the int value to be read
    @return boolean value indicating wether a value has been read or
    dnot. 
*/ 
int readval(int &rflag,char *line,int &val);

#if 0 // This is already defined in 'util2.h'
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Reads a text hash table from a filestack. 
    @author M. Storti
    @param fstack (input) The filestack from which the hash table is
    read. 
    @param (output) The hash table that has been read. 
*/ 
int read_hash_table(FileStack *& fstack,TextHashTable *& thash);
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Prints a vector with ``rotary save'' mechanism. 
    @author M. Storti
    @param filenamepat (input) The pattern that generates the
    filenames to which the vectors are written. Must contain #%d#
    @param x (input) The vector to be written
    @param dofmap (input) The dofmap of the problem. (Gives values for
    fixed node/field combinations - Dirichlet boundary conditions.)
    @param time_data (input) The corresponding time instant. (Gives values for
    fixed node/field combinations - Time dependent Dirichlet boundary
    conditions .)
    @param j (input) The actual time step
    @param nrec (input) The number of records allowed in each file. 
    @param nfile (input) The number of files allowed. 
*/ 
void print_vector_rota(const char *filenamepat,const Vec x,const
		       Dofmap *dofmap,const TimeData *time_data,
		       const int j,const int nsave,const int nrec,
		       const int nfile);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Parses a line of props of the form 'prop1[len1] prop2[len2] ...'
    Some properties may not have a length in which case we asume 1.
    @author M. Storti
    @param line (input) the props line
    @param prop_name (output) a vector of strings containing the
    property names
    @param prop_len (output) a vector of ints containing the length of
    the property 
*/ 
void parse_props_line(const char *line,vector<string> &prop_name,vector<int>
		      &prop_len);

#endif
