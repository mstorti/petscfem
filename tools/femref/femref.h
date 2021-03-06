// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: femref.h,v 1.63.70.1 2007/03/25 03:19:40 mstorti Exp $
#ifndef PETSCFEM_FEMREF_H
#define PETSCFEM_FEMREF_H

#include <list>
#include <map>
#include <deque>
#define Mesh __Mesh__
#include <src/dvector.h>
#include <src/dvector2.h>
#undef Mesh
#include <src/generror.h>
#include <src/linkgraph.h>
#include "./tree.h"

// using namespace __gnu_cxx;

extern MD5Hasher hasher;

/// This is a value that a real node can't take. 
#define NULL_NODE INT_MAX

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Represents a geometric object `outside' its container.
    Basically, it represents a sequence of nodes and a 
    topological shape (template). */ 
class GeomObject {
public:
  /** Types of objects may be identified by a #Type#
      value, or as a pointer to a #Template# class. 
   */ 
  enum Type { NULL_TYPE=0, 
	      EdgeRefNodeT, 
	      OrientedEdgeT, EdgeT, 
	      OrientedTriT, TriT,
	      OrientedTetraT, TetraT};
  /// Contains information about a shape
  class Template {
  public:
    // friend class GeomObject;
    // friend class Mesh;
    /// Contains the permutations that leave the shape invariant
    dvector<int> perms_v;
    /// Number of nodes, dimension, number of permutations
    int size_m,dim_m,nperms_m;
    /// The corresponding enum #Type#. 
    GeomObject::Type type;
    /// The name of the shape
    const char *label;
    virtual ~Template() {}
    /// Ctor from info
    Template(int sz,GeomObject::Type t,
	     int dim_a,int nperms_a,
	     const int *perms,const char *label_a);
    /// The number of objects of type #t# this shape has. 
    virtual int size(Type t) const { return 0; }
    /// Local nodes connected to subobject #j# of type #t#. 
    virtual const int* nodes(Type t,int j) const { assert(0); return NULL; }
  };
  /// Ptr to shape structure
  Template *go_template;
  /** This function returns the shape ptr for a given
      enum #Type# */ 
  static Template* get_template(Type t);
private:
  /// The check-sum for this object, has been converted to canonical?
  int cs, canonical, hash_value;
  /// The nodes for this object
  dvector<int> nodes_m;
public:
  /// Number of nodes for this object
  int size() const { return go_template->size_m; }
  /// The number of objects of type #t# this shape has. 
  int size(Type t) const { return go_template->size(t); }
  /// The checksum for this object
  int csum() const { return cs; }
  /// The hash value for this object
  int hash_val() const;
  /// Dimension of the object
  int dim() const { return go_template->dim_m; }
  /// Number of permutations that leave this object invariant. 
  int nperms() const { return go_template->nperms_m; }
  /// Reset to the empty object
  void clear() { go_template=NULL; nodes_m.clear(); canonical=0; cs=0; }
  /// Set to an object with the specific info
  void init(Type t,const int *nodes_a=NULL);
  /// Set to an object with #local_nodes#
  /// from array #global_nodes#
  void init(Type t,const int *local_nodes,
	    const int *global_nodes);
  /// Local nodes for #j# permutation 
  const int* perm(int j) const {
    return &(go_template->perms_v.e(j,0)); 
  }
  /// Returns the #enum# type. 
  Type type() const { return go_template->type;}
  /// Constructor of empty object
  GeomObject() : go_template(NULL), canonical(0) { }
  /// Constructor with info
  GeomObject(Type t,const int *nodes_a=NULL);
  /** Bring to canonical ordering by permuting the
      nodes.  The canonical ordering is that
      permutation of the nodes that gives the lowest
      lexicographical order of the sequence of nodes.
      The sequence of nodes in canonical order is
      a unique identification of the object. 
   */ 
  void make_canonical();
  /// Compare two objects
  bool equal(GeomObject &go);
  /// Finds the subobject that matches object #go#. 
  int find(GeomObject &go) const;
  /// Prints the object
  void print(const char*s = NULL) const;
  /// The nodes of this object
  const int* nodes() const { return nodes_m.buff(); }
  /** Set #go# to be the #j#-th subobject of type #t# in
      #*this# */ 
  void set(Type t,int j,GeomObject &go) const;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class Splitter {
public:
  /// Number of subobjetcs of type #t# in this splitting.
  virtual int size(GeomObject::Type t) const { assert(0); return -1; }
  /// Number of refined nodes this splitting creates
  int nref_nodes() const { assert(0); return -1; }
  /** Refined nodes #indx# is created from #nnod# nodes 
      passed by #nodes#. 
      @param indx (input) the refined node index. 
      @param tmpl (output) the type of the refined 
      node (may be edge, face, ...)
      @param nnod (output) the number of nodes #indx# 
      refined node depends. 
      @param nodes (output) the nodes that spawn #indx#. */ 
  void ref_node(int indx,
		const GeomObject::Template *&tmpl,
		int &nnod, const int *&nodes) 
    const { assert(0); }
  /// Local nodes connected to subobject #j# of type #t#. 
  virtual const int* 
  nodes(GeomObject::Type t,int j) const { assert(0); return NULL; }
  /// Total number of subobjetcs
  virtual int size() const { assert(0); return -1; }
  /// Local nodes connected to subobject #j# and type
  virtual const int 
  *nodes(int j,GeomObject::Type &t) const { assert(0); return NULL; }
  virtual ~Splitter() {}
};

typedef double
(*RefineFunction)(GeomObject &go,const double *xnod);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** The data to be stored for a node. Data for
    derived nodes (created during refinement) are
    computed through recursive use of a 
    #NodeInfoCombineFunction#. */ 
class NodeInfo {
public:
  virtual ~NodeInfo()=0;
  virtual void print() { };
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** The type for the container that
    stores NodeInfo objects */ 
typedef 
map<int,NodeInfo *> NodeInfoMapT;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Model class that combines the data
    for two nodes and then creates the data for
    the combined node. */ 
class NodeCombiner {
public:
  virtual ~NodeCombiner()=0;
  virtual void 
  combine(int tag,int n,
	  const int *nodes,
	  int newnode,
	  NodeInfoMapT &node_info_map)=0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** A mesh is a container of geometrical objects
    linked (as in a graph).  Two objects are linked
    if they share a node. For reasons of efficiency,
    meshes are not stored by storing the whole graph. */ 
class Mesh {
public:
  /** An iterator is a position in the
      graph. Basically it is the number of element
      in the mesh plus the number that identifies
      the subobject inside the element. Iterators do
      not represent uniquely the object,
      i.e. several iterators can point to the same
      object. */ 
  class iterator {
  public:
    // friend class Mesh;
    /// The big object number in the mesh
    int obj;
    /** The type of the geometric object
	identified by this iterator */ 
    GeomObject::Type t;
    /// The position of the subobject in the large object
    int subobj;
    /// Ctor from the data
    iterator(int obj_a,GeomObject::Type ta,int subobj_a)
      : obj(obj_a), t(ta), subobj(subobj_a) { }
    /// Ctor of empty iterator
    iterator() : obj(-1), t(GeomObject::NULL_TYPE), 
		 subobj(-1) { }
    /// Set to data 
    void set(int obj_a,GeomObject::Type ta,int subobj_a) {
      obj = obj_a; t = ta; subobj = subobj_a; 
    }
  };
  virtual ~Mesh() {}
  /** Set object #go# to the object pointed
      by iterator #it# */
  virtual void set(iterator it,GeomObject &go)=0;
  /** Returns one iterator that points to the 
      specified object. Remember that may be many. */ 
  virtual iterator find(GeomObject &go)=0;
  /** Finds the list of all iterators that point to
      object #go#. */ 
  virtual void find(GeomObject &go,list<iterator> &its)=0;
  /** Get list #adj# of iterators of shape #t# that
      that are connected to iterator #t#  */ 
  virtual void get_adjacency(iterator it,GeomObject::Type t,
			     list<iterator> &adj)=0;
  /// Determine whether this iterator is empty. 
  virtual bool is_end(iterator it)=0;
};

/** This is a container of objects formed by 
    all large objects of shape #tmpl#. */ 
class UniformMesh : public Mesh {
private:
  /// The coordinates of the nodes
  dvector<double> coords;
  /** The connectivities of the mesh. List of
      nodes of each large object. */
  dvector<int> connec;
  /** Adjacency list node-to-element. Stores the list of
      elements connected to a given node. */ 
  dvector<int> n2e;
  /** Pointers to ranges in #n2e#, elements connected 
      to node #j# are in range #[n2e_ptr[j],n2e_ptr[j+1])#. 
      (#n2e_ptr# is size #nnod+1#). */
  dvector<int> n2e_ptr;
  /// The dimension of the space
  int ndim;
  /** Number of nodes, number of elements, number of
      nodes connected to an element */
  int nnod, nelem, nel;
  /// The ptr to the shape object
  const GeomObject::Template *tmpl;
  struct ElemRefNode {
    const Splitter *splitter;
    int so_indx;
  };
  typedef tree<ElemRefNode> ElemRef;
  /// The refinement of each element
  dvector<ElemRef *> elem_ref;

  /// An iterator to identify refined nodes
  class RefNodeIterator {
  public:
    ElemRef::cell *c;
    int indx;
  };

  /** Auxiliary function that searches the iterator
      for a given object. If #its==NULL# then
      returns the first iterator found in #it#. If
      #its!=NULL# then returns the list of iterators
      in #its#.
      @param go (input) the geometric object for which 
      to return the list of iterators. 
      @param its (output) if not null, then return the list 
      of iterators in its. 
      @param it (output) if #its=NULL# then return the 
      first iterator found here. */ 
  void find(GeomObject &go,list<iterator> *its,iterator &it);
  typedef multimap<int,RefNodeIterator> hash2it_t;
  /// For a given `NodeRef' hash value, gives an iterator for it
  hash2it_t hash2it;
  /// Stores the correspondence between hash values and nodes
  map<int,int> hash2node;
  /** Nodes obtained by refinement are in the range 
      #[nnod,last_ref_node)# */ 
  int last_ref_node;
  /** Combines two base nodes in order to get a new
      node number. */ 
  int combine(int n1,int n2);

public:
  /// Ctor from dimensions and shape
  UniformMesh(const GeomObject::Template &tmpl_a,int ndim_a); 
  /// Dtor
  ~UniformMesh(); 
  const GeomObject::Template *tmplt() { return tmpl; }
  /// Init mesh from connectivities
  void set_conn(const dvector<int> &icone,
	    int base=0);
  /// Read the mesh from specified files
  void read(const char *node_file,
	    const char *conn_file,
	    int base=0);
  /** Sets #go# to object pointed by #it#. 
      @param it (input) iterator to geometric object in mesh
      @param go (output) the opject pointed by #it# */ 
  void set(iterator it,GeomObject &go);

  /** Builds subobject #sgo# to be the #indx#-th subobject
      of #go# when splitted with splitter #s#. 
      @param elem (input) the element index
      @param go (output) the object to be constructed
      @param ref_nodes (output) the nodes for #sgo# that
      @param comb (input) the function used to 
      combine NodeInfo objects
      @param node_info_map (input/output) adds new #NodeInfo#
      objects to this container */ 
  void set(int elem,
	   GeomObject &sgo,
	   list<int> &ref_nodes,
	   NodeCombiner *node_comb,
	   NodeInfoMapT *node_info_map);  

  /** Builds subobject #sgo# to be the #indx#-th subobject
      of #go# when splitted with splitter #s#. 
      @param go (input) the "father" object
      @param s (input) the splitter for the father object
      @param indx (input) the subobject index for this object
      @param sgo (output) the object to be constructed
      @param ref_nodes (output) the nodes for #sgo# that
      @param comb (input) the function used to 
      combine NodeInfo objects
      @param node_info_map (input/output) adds new #NodeInfo#
      objects to this container */ 
  void set(const GeomObject &go,
	   const Splitter *s,
	   int indx,
	   GeomObject &sgo,
	   list<int> &ref_nodes,
	   NodeCombiner *node_comb,
	   NodeInfoMapT *node_info_map);  
  /** Finds the _first_ iterator to object #go#. 
      If there isn't, returns the empty iterator. 
      @param go (input) the object to find
      return first iterator that points to #go# */ 
  iterator find(GeomObject &go);
  /** Finds _all_ the iterators that point to object #go#. 
      @param go (input) the object to find
      @param its (output) list of iterators that
      point to #go# */ 
  void find(GeomObject &go,list<iterator> &its);
  /** Get the list of iterators of type #t# 
      adjacent to object pointed by #it#. 
      @param it (input) iterator to object for which to get 
      the adjacency list
      @param t (input) find adjacent objects of type #t#
      @param adj (output) the list of adjacent iterators */ 
  void get_adjacency(iterator it,GeomObject::Type t,
		     list<iterator> &adj) { 
    assert(0); // not implemented yet!
  }
  /** Predicate for determine if an iterator is the empty one. 
      @param it (input) iterator to check
      @return #true# if #it# is empty, #false# otherwise. */ 
  bool is_end(iterator it) { 
    return it.obj<0 || it.t==GeomObject::NULL_TYPE 
      || it.subobj<0; 
  }
  /** Refines mesh acording to #RefineFunction#
      @param rf (input) function that indicates 
      the desired mehs size (#h#) for this element. */ 
  void refine(RefineFunction rf);

  /** Given an iterator to a #RefNode# constructs in 
      #go# the corresponding geometrical object. 
      @param it (input) the iterator
      @param go (input) the geometrical object (ogf type 
      refined node). */ 
  void set(RefNodeIterator it,
	   GeomObject &go);

  /** Adds new refined node to the hash
      with given arguments */
  void add_refined_node(ElemRef::iterator q,
			int j,int node_hash);

private:
  void hash_insert(int k,RefNodeIterator q);

public:
  /** Each element in the refinement stack is
      of this type. */ 
  struct RefPathNode {
    /// The geometric object
    GeomObject go;
    /** The splitter for this object. Actually it
	is a node in the refinement tree for this
	element */ 
    ElemRef::iterator splitter;
    /** The sibling index in the sibling list for
	this object */
    int so_indx;
    /** List of nodes that have been created
	at this level */ 
    list<int> ref_nodes;
  };

  friend class visitor;

  enum VisitMode { Natural=0, BreadthFirst };
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Refines mesh acording to #RefineFunction#
      @param rf (input) function that indicates 
      the desired mehs size (#h#) for this element. */ 
  class visitor {
  private:
    /// The mesh we are visiting
    UniformMesh *mesh;
    /// The internal tree for each base element
    ElemRef *etree_p;    
    /// The element we are currently visiting
    int elem;
    /// Pop the deepest level object in the refinement stack
    void pop();
    /// Go to next basic element. This is a callback
    /// for defining breadth first and natural order traversing
    /// algorithms. 
    void next_elem();
    bool end_elem();
    int nvisited;
    dvector<char> visited;
    deque<int> element_stack;
    void pop_elem();
  public: 
    VisitMode visit_mode;
    /// The type for the refinement stack
    typedef list<RefPathNode> RefStackT;
    /// Flags whether we print the elements as they are visited
    int trace;
    /// Functions used to combine #NodeInfo# at nodes
    NodeCombiner *node_comb;
    /// Stores `NodeInfo' objects
    NodeInfoMapT node_info_map;
    /// Ctor.
    visitor();
    /// Stack containing the elements in the refinement tree
    RefStackT ref_stack;
    /// Inits the visitor to the first element of #mesh#
    void init(UniformMesh &mesh_a,int elem=0);
    /// Inits the visitor to the base element of index #elem#
    void init(int elem);
#if 0
    /** Inits the visitor to the base element of the
	first element. */
    void init(UniformMesh &mesh_a);
#endif
    /** Pass to the following subobject of the 
	same element. Return false if reached the end. */
    void so_next();
    /// Have we passed the last subobject of this element?
    bool so_end();
    /** Pass to the following subobject of the 
	mesh. Return false if reached the end of the mesh. */
    void next();
    /** Pass to the following subobject of the 
	mesh, down to level #level#. 
	Return false if reached the end of the mesh. */
    void next(int level);
    /// Have we passed the last subobject of the mesh?
    bool end();
    /** Pass to the following subobject of the 
	mesh, at this level or higher. */
    bool level_next();
    /** Pass to the following subobject of the 
	element, at this level or higher. */
    bool level_so_next();
    /** Is this node a leave in the refinement tree? */ 
    bool is_leave();
    /** Refine this element according to #s# */ 
    void refine(const Splitter* s);
    /** The refinement level for this node. */ 
    int ref_level();
    /** Return the currently visited element. */
    int elem_indx() const { return elem; }
  };
  friend class LinearCombiner;
};


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class NodeInfoRef : public NodeInfo {
public:
  vector<double> coords;
  ~NodeInfoRef() {  }
  void print();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class LinearCombiner : public NodeCombiner {
public:
  UniformMesh *mesh;
  ~LinearCombiner() { } 
  void 
  combine(int tag,
	  int n,const int *nodes,
	  int new_node,
	  NodeInfoMapT &node_info_map);
};

extern LinearCombiner linear_combiner;

struct GetSurfCtx {
  UniformMesh *mesh;
  int ndim, nnod;
  GetSurfCtx() : mesh(NULL) { }
  ~GetSurfCtx() { if (mesh) delete mesh; }
};

void getsurf(GetSurfCtx &ctx,
	     const dvector<int> &icone,
	     dvector<int> &surf_con,
	     dvector<int> &surf_nodes, 
	     int base, int verbose);

void comp_matrices(GetSurfCtx &ctx,
		   const dvector<int> &surf_con,
		   const dvector<int> &surf_nodes, 
		   dvector<double> &x, 
		   dvector<double> &surf_mass,
		   dvector<double> &node_mass,
		   int verbose);

void 
elem2nod_proj(GetSurfCtx &ctx,
	      const dvector<int> &icone,
	      const dvector<double> &elem_mass,
	      const dvector<double> &node_mass,
	      const dvector<double> &ue,
	      dvector<double> &un);

void
nod2elem_proj(GetSurfCtx &ctx,
	      const dvector<int> &icone,
	      const dvector<double> &un,
	      dvector<double> &ue);

#if 0
void 
fem_smooth(GetSurfCtx &ctx,
	   const dvector<int> &surf_con,
	   const dvector<double> &surf_mass,
	   const dvector<double> &node_mass,
	   const dvector<double> &u,
	   dvector<double> &us,
	   int niter,
	   int verbose);
#endif

#endif
