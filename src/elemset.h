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

#ifndef ELEMSET_H
#define ELEMSET_H

#include "arglist.h"
#include "libretto.h"
#include <glib.h>

enum ElemsetIteratorMode {
  ALL                  = 0x00001,
  INCLUDE_GHOST        = 0x00010,
  DO_NOT_INCLUDE_GHOST = 0x00100
};

class ElementList;
class ElementIterator;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Stores sets of elements with similar types and properties. 
    @author M. Storti
*/ 
class Elemset {
public:
  // type of element
  char *type;
  /// table of connectivities
  int *icone;
  ///number of elements in the elemset
  int nelem;
  ///number of nodes per element
  int nel;
  ///number of degrees of freedom per node
  int ndof;
  /// mesh partitioning as computed by Metis 
  int *epart;
  /// flag indicating whether this  is the fat elemset or not
  int isfat;
  /// number of elements in this processor
  int nelem_here;
  /// number of double properties in the per-element properties table
  int nelprops;
  /// number of integer properties in the per-element properties table
  int neliprops;
  /// number of additional double properties 
  int nelprops_add;
  /// number of additional integer properties
  int neliprops_add;
  /// double per-element properties table
  double *elemprops;
  /// int per-element properties table
  int *elemiprops;
  /// additional double per-element properties table
  double *elemprops_add;
  /// additional int per-element properties table
  int *elemiprops_add;
  /// This is a ``loccker'' for each element
  void **local_store;
  /// properties hash table 
  TextHashTable *thash;
  /// hash of properties in the per element double properties table
  GHashTable *elem_prop_names;
  /// hash of properties in the per element int properties table
  GHashTable *elem_iprop_names;

  /** Assembles residuals, matrices, scalars and other quantities. 
      @author M. Storti
      @param retval (output) here returns contributions vectors
      (residuals), and matrices
      @param nodedata (input) vector with properties per node
      @param locst (input) state vectors localized to elements 
      @param locst2 (input) same for alternative  state vector (this may be
      used in temporal integration), etc...
      @param dofmap (input) maps the node/field representation to the vector
      state
      @param ijob (input) tells the global assemble function  what kind of job
      should be done (assemble vector, profile or matrix). To be
      processed by the global assemble function. 
      @param jobinfo (input) tells the routine element what kind of
      matrix or vector has to be assembled
      @param myrank (input) identifies this processor
      @param el_start (input) low end of the range of elements to be
      processed
      @param el_last (input) high end of the range of elements to be
      processed
      @param iter_mode (input) include or not ghost elements
  */
  virtual int assemble(arg_data_list &arg_datav,Nodedata *nodedata,Dofmap *dofmap,
		       char *jobinfo,int myrank,
		       int el_start,int el_last,int iter_mode,
		       const TimeData *time_data) {
    
    printf("assemble: not known Elemset\n"); exit(1);};

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Ask the elemset if it should be processed for this jobinfo. 
      @author M. Storti
      @param jobinfo (input) The name of the task.
      @param answer (output) 1 = process (0 = do not process) this
      elemset. 
  */ 
  virtual int ask(char *jobinfo,int &skip_elemset) {

    // By default process the elemset. 
    skip_elemset = 0;
    return 0;
  };
    

  /// dynamic aray with ghost elements
  Darray *ghost_elems;
  /// print info of this elemset
  void print();

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Builds localized vectors (in node/field) representation of state
      vectors. 
      
      @author M. Storti
      @param nel (input) number of nodes per element
      @param ndof (input) number of degrees of freedom per node
      @param dofmap (input) the dofmap of the mesh
      @param locst (output) localized vector to be assembled
      @param myrank (input) identifies this processor
      @param el_start (input) low end of the range of elements to be
      processed
      @param el_last (input) high end of the range of elements to be
      processed
      @param iter_mode (input) include or not ghost elements
      @param time_data (input, def=NULL) an external parameter in order to compute
      external boundary conditions, etc...
  */ 
  int download_vector(int nel,int ndof,Dofmap *dofmap,
		      arg_data &argd,
		      int myrank,int el_start,int el_last,int iter_mode,
		      const TimeData *time_data=NULL);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Put localized values (nod/field representation) returned by
      elemset assembles in the global residual vector. 
      
      @author M. Storti
      @param nel (input) number of nodes per element
      @param ndof (input) number of degrees of freedom per node
      @param dofmap (input) the dofmap of the mesh
      @param retval (input) localized values (nod/field
      representation) returned by elemset assembles 
      @param A (input/output) assemble matrix contributions on it (if
      ijob/jobinfo indicates so). 
      @param da (output) assemble matrix profile on this (if
      ijob/jobinfo indicates so). 
      @param vec (output) assemble vector contributions on it (if
      ijob/jobinfo indicates so).
      @param myrank (input) identifies this processor
      @param el_start (input) low end of the range of elements to be
      processed
      @param el_last (input) high end of the range of elements to be
      processed
      @param iter_mode (input) include or not ghost elements
      @param klocc (input) if computing difference finite jacobian, is
      the local node perturbed 
      @param klocc (input) if computing difference finite jacobian, is
      the local d.o.f. perturbed 
  */
  int upload_vector(int nel,int ndof,Dofmap *dofmap,
		    int options,arg_data &argd,int myrank,
		    int el_start,int el_last,int iter_mode,
		    int klocc=0,int kdofc=0);

  /** Return the ``locker'' for the element, that is a (void *) where tu
      put data. Only for elements local to this processor. 
      @author M. Storti
      @param local_elem (input) the index 
      (local to the processor) of the element
      @return the address (void *) of the locker
  */
  void *& local_store_address(int local_elem) {
    return local_store[local_elem];
  }

  friend class ElementList;

};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Global assemble function. Loops over all elemsets and assembles
    element contributions to amtrices, vectors and profiles. 
    Same as assemble, but for computine profiles.
    @author M. Storti
    @param mesh (input) the mesh to be processed
    @param argl (input/output) list of arguments, (vector matrices) on input
    and output. 
    @param dofmap (input) the dofmap of the mesh
    @param jobinfo (input) tells the routine element what kind of
    matrix or vector has to be assembled
    @param time_data (input, def=NULL) an external parameter in order to compute
    external boundary conditions, etc...
    @return error code
*/
int assemble(Mesh *mesh,arg_list argl,Dofmap *dofmap,char *jobinfo,
	     const TimeData *time_data=NULL);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Const iterator for looping over elements in an elemset.
class ElementList {
private:
  /// Pointer to the elemset
  const Elemset *elemset;
  /// Start range of values while iterating
  int first;
  /// Ends range of values while iterating
  int last;
  /// Iteration mode
  ElemsetIteratorMode mode;
public:
  /// Constructor
  ElementList(const Elemset *elemset_,
	      int first_,int last_,ElemsetIteratorMode mode_) :
    elemset(elemset_), first(first_), last(last_), mode(mode_) {};
  /// Allow iterator to access this class
  friend class ElementIterator;
  /// Begin of chunk
  ElementIterator begin(void) const;
  /// End of chunk
  ElementIterator end(void) const;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Iterators for looping over elements in an elemset.
class ElementIterator {
private:
  /// The list to which it belongs
  const ElementList *elemlist;
  /// Rank of element in the chunk
  int rank_in_chunk;
  /// Rank of element in the elemset
  int rank_in_elemset;
public:
  /// Constructor
  ElementIterator(const ElementList *el, const int rie,const int ric) 
    : elemlist(el), rank_in_elemset(rie), rank_in_chunk(ric) {};
  /// Prefix increment operator
  ElementIterator & operator++(void);
  /// Postfix Increment operator
  ElementIterator operator++(int);
  /// Not equal operator
  int operator!=(const ElementIterator &other) const;
  /// Return position in chunk
  void position(int &pos_in_elemset, int &pos_in_chunk) {
    pos_in_elemset = rank_in_elemset; pos_in_chunk = rank_in_chunk;
  }
  /// Returns true/false whether the current element is
  /// in the restricted list or not
  int is_valid() const;
  /// Advances iterator to 
  /// in the restricted list or not
  void advance_to_valid();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** @name compute\_prof package */
//@{
/** nodes to store the profile of a matrix. 
    @author M. Storti
*/ 
class Node {
public:
  /// points to next dof connected to this one
  int next,
    /// the connected dof
    val;
  /// constructor
  Node(int next_=-1 ,int val_=0 ) : next(next_), val(val_) {};
  /// printer
  void print(void) {printf(" (%d,%d)\n",next,val);};
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Inserts a node in the profile. 
    @author M. Storti
    @param da (input/output) dynamic array containing the profile
    @param j (input) row index
    @param k (output) column index
*/ 
void node_insert(Darray *da,int j,int k);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Compare function for sort and then binary search on the profile. 
    @author M. Storti
    @param left (input) pointer to left element to be compared
    @param right (input) pointer to right element to be compared
    @param args (not used) as required by libretto routines
*/ 
int int_cmp (const void *left,const void *right, void *args);
//@}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Flags if the current element is to be computed his contribution. 
    Basically it process all elements in the current
    processor. However, depending on iter\_mode it may include also
    ghot-elements. 

    @author M. Storti
    @param iele (input) element to be tested
    @param elemset (input) this elemset 
    @param myrank (input) this processor index
    @param iter\_mode (input) flags whether or not include
    ghost-elements. 
    @return boolean flag indicating whether process or not this element
*/ 
int compute_this_elem(const int & iele,const Elemset *elemset,const int & myrank,
		      int iter_mode);

#define ELEMSET_CLASS(name) \
class name : public Elemset { \
public: \
  int assemble(arg_data_list &arg_data_v,Nodedata *nodedata, \
	       Dofmap *dofmap,char *jobinfo,int myrank, \
	       int el_start,int el_last,int iter_mode); \
}

#define ASSEMBLE_FUNCTION \
  int assemble(arg_data_list &arg_data_v,Nodedata *nodedata, \
	       Dofmap *dofmap,char *jobinfo,int myrank, \
	       int el_start,int el_last,int iter_mode, \
	       const TimeData *)

#define ASK_FUNCTION int ask(char *jobinfo,int &answer)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Loops many times executing assmble() with the same arguments and
    issuing a performance indication. 
    @author M. Storti
    @param mesh (input) the mesh to be processed
    @param argl (input/output) list of arguments, (vector matrices) on input
    and output. 
    @param dofmap (input) the dofmap of the mesh
    @param jobinfo (input) tells the routine element what kind of
    matrix or vector has to be assembled
    @param time_data (input, def=NULL) an external parameter in order to compute
    external boundary conditions, etc...
    @return error code
*/ 
int measure_performance_fun(Mesh *mesh,arg_list argl,
			Dofmap *dofmap,char *jobinfo,const TimeData
			*time_data=NULL);

typedef void 
NewAssembleFunction(arg_data_list &arg_datav,const Nodedata *nodedata,const Dofmap *dofmap,
		    const char *jobinfo,const ElementList &elemlist,
		    const TimeData *time_data);

/** This is an adaptor to the old Elemset class
    @author M. Storti
*/
class NewElemset : public Elemset {
  /// This is the adaptor to the old assemble function.
  int assemble(arg_data_list &arg_datav,Nodedata *nodedata,Dofmap *dofmap,
	       char *jobinfo,int myrank,
	       int el_start,int el_last,int iter_mode,
	       const TimeData *time_data);
  /// The new assemble function
  virtual void 
  new_assemble(arg_data_list &arg_datav,const Nodedata *nodedata,const Dofmap *dofmap,
	       const char *jobinfo,const ElementList &elemlist,
	       const TimeData *time_data) {
    printf("assemble: not known New Elemset\n"); exit(1);
  };
};

#if 0
#define NEW_ASSEMBLE_FUNCTION \
  int assemble(arg_data_list &arg_data_v,Nodedata *nodedata, \
	       Dofmap *dofmap,char *jobinfo,int myrank, \
	       int el_start,int el_last,int iter_mode, \
	       const TimeData *)
#endif

#endif
