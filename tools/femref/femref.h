// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: femref.h,v 1.1 2004/11/16 02:59:56 mstorti Exp $
#ifndef PETSCFEM_FEMREF_H
#define PETSCFEM_FEMREF_H

#include <src/dvector.h>
#include <src/dvector2.h>


#if 0
class ObjectTable

class NodeSet {
 public:
  /// Ctor 
  NodeSet(int ndim_a) : ndim(ndim_a), max(0) { }
  /** A node is identified by an index in `node_list'. 
      All indices between in range #[0,max-1)# that are not in 
      #invalid#, are nodes. **/
  class NodeHandle : public int { };
  /// Adds a new node with those specific coordinates
  NodeHandle add(double *coords);
  /// Remove that nodes
  void remove(NodeHandle node);
  
 private:
  /// The dimension of the space
  int ndim;
  /* Stores the coordinates.  A node is identified by an index in
     `node_list'.  All indices in range #[0,max-1)# that are not in
     #invalid#, are nodes. **/
  dvector<double> coords;
  /// Max node added so far
  int max;
  /// These are `holes' in the container
  set<int> invalid;
};
#endif

class Mesh {
 public:
  class Object {};
  
 private:
  dvector<double> xnod();
  dvector<int> icone();
};

#endif
