//__INSERT_LICENSE__
// $Id: visitor.cpp,v 1.3 2004/12/19 14:28:49 mstorti Exp $

#include <string>
#include <list>
#include "./hasher.h"

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
  list<RefPathNode>::iterator ws, 
    w = ref_stack.begin();
  // here visit w...
  ElemRef::iterator q, qs, qrsib, qfather;
  q = w->splitter;
  int j = w->so_indx;
  printf("level %d, sibling %d,",ref_stack.size()-1,j);
  w->go.print("");

  if (q != etree_p->end()) {
    // `q' is a regular node for splitters (sure
    // it isn't a leave for GO's). Follow to the
    // left child (int the GO tree). 
    qs = q.lchild();
    ref_stack.push_front(RefPathNode());
    ws = ref_stack.begin();
    w = ws; w++;
    // Build `ws' from GeomObject `w' (parent) and
    // splitter `s' and subobject index `j'
    mesh->set(w->go,q->splitter,j,ws->go);
    ws->go.make_canonical();
    ws->so_indx = 0;
  } else {
    // `q' is a leave for GO's (sure it isn't a
    // regular node for splitters). Try to find a
    // right sibling, or a father that has a right
    // sibling
    while (true) {
      // Check if we are at the root
      if (ref_stack.size()<=1) 
	return false;
      
      w = ref_stack.begin(); w++;
      qfather = w->splitter;
      assert(qfather != etree_p->end());
      const Splitter *s = qfather->splitter;

      ref_stack.pop_front();

      if (j<s->size()-1) {
	int jsib = j+1;
	ref_stack.push_front(RefPathNode());
	ws = ref_stack.begin();
	w = ws; w++;
	// Build `ws' from GeomObject `w' (parent) and splitter `s'
	// and subobject index `jsib'
	mesh->set(w->go,s,jsib,ws->go);
	ws->go.make_canonical();

	// Find next node on the splitting tree or end()
	qrsib = ws->splitter;
	while (qrsib != etree_p->end()) {
	  if (qrsib->so_indx >= jsib) break;
	  qrsib++;
	}

	// Push new state in the stacks
	ws->splitter = qrsib;
	ws->so_indx = jsib;
	return true;
      }
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void UniformMesh::visitor::
init() { assert(0); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool UniformMesh::visitor::
next() { assert(0); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool UniformMesh::visitor::
so_end() { assert(0); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool UniformMesh::visitor::
end() { assert(0); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void UniformMesh::visitor::
refine() {
  assert(is_leave());
  list<RefPathNode>::iterator 
    w = ref_stack.begin();
  ElemRef::iterator 
    q = w->splitter;
  q = etree_p->insert(q,ElemRefNode());
  w->splitter = q;
  int j = w->so_indx;
  // The splitter should be returned
  // by the refinement function
  q->splitter = &Tetra2TetraSplitter;
  q->so_indx = j;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool UniformMesh::visitor::
is_leave() {
  list<RefPathNode>::iterator 
    w = ref_stack.begin();
  return w->splitter == etree_p->end();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int UniformMesh::visitor::
ref_level() { return ref_stack.size()-1; }
