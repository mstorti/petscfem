//__INSERT_LICENSE__
// $Id: visitor.cpp,v 1.1 2004/12/18 22:34:01 mstorti Exp $

#include <list>

using namespace std;

#include "./femref.h"
#include "./gtemplates.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void UniformMesh::visitor::init(UniformMesh &mesh_a,int elem) {
  ref_stack.clear();
  mesh  = &mesh_a;
  etree_p = mesh->elem_ref.e(elem);
  ElemRef::iterator q = etree_p->begin();
  
  ref_stack.push_front(RefPathNode());
  RefPathNode &top = ref_stack.front();
  top.go.init(mesh->tmpl->type,&mesh->connec.e(elem,0));
  top.go.make_canonical();
  top.splitter = etree_p->begin();
  top.so_indx = 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool UniformMesh::visitor::so_next() {
  RefPathNode &w = ref_stack.begin();
  // here visit w...
  ElemRef::iterator q, qs, qrsib, qfather;
  q = split_stack.front();
      int j = split_indx_stack.front();
      printf("level %d, sibling %d,",go_stack.size()-1,j);
      w->print("");

      if (q==etree.end()) {
	// May be refine this element?
#define REF_LEVEL 2
	  // Here refine (eventually) by inserting a
	  // child in the position `q'.
	  // Refinement criterion.
	  bool refine_this = go_stack.size()<=REF_LEVEL;
	  if (refine_this) {
	    q = etree.insert(q,ElemRefNode());
	    // The splitter should be returned
	    // by the refinement function
	    q->splitter = &Tetra2TetraSplitter;
	    q->so_indx = j;
	    split_stack.front() = q;
	  }
      }

      if (q != etree.end()) {
	// `q' is a regular node for splitters (sure
	// it isn't a leave for GO's). Follow to the
	// left child (int the GO tree). 
	s = q->splitter;
	qs = q.lchild();
	split_stack.push_front(qs);
	go_stack.push_front(GeomObject());
	ws = go_stack.begin();
	// Build `ws' from GeomObject `w' (parent) and
	// splitter `s' and subobject index `j'
	set(*w,s,j,*ws);
	ws->make_canonical();
	split_indx_stack.push_front(0);
      } else {
	// `q' is a leave for GO's (sure it isn't a
	// regular node for splitters). Try to find a
	// right sibling, or a father that has a right
	// sibling
	while (true) {
	  // Check if we are at the root
	  if (split_stack.size()<=1) {
	    done = true; break;
	  }
	  list<ElemRef::iterator>::iterator 
	    qit = split_stack.begin();
	  qfather = *(++qit);
	  assert(qfather != etree.end());
	  s = qfather->splitter;
	  j = split_indx_stack.front();

	  go_stack.pop_front();
	  split_stack.pop_front();
	  split_indx_stack.pop_front();
	  if (j<s->size()-1) {
	    int jsib = j+1;
	    go_stack.push_front(GeomObject());
	    ws = go_stack.begin();
	    w = ws; w++;
	    // Build `ws' from GeomObject `w' (parent) and splitter `s'
	    // and subobject index `jsib'
	    set(*w,s,jsib,*ws);
	    ws->make_canonical();

	    // Find next node on the splitting tree or end()
	    qrsib = q;
	    while (qrsib != etree.end()) {
	      if (qrsib->so_indx >= jsib) break;
	      qrsib++;
	    }

	    // Push new state in the stacks
	    split_stack.push_front(qrsib);
	    split_indx_stack.push_front(jsib);
	    break;
	  }
	}
      }

}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void UniformMesh::visitor::init() { assert(0); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool UniformMesh::visitor::next() { assert(0); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool UniformMesh::visitor::so_end() {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool UniformMesh::visitor::end() { assert(0); }
