//__INSERT_LICENSE__
// $Id: visitor.cpp,v 1.2 2004/12/18 22:45:26 mstorti Exp $

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
  list<RefPathNode>::iterator ws, w 
    = ref_stack.begin();
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
#if 0
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
#endif
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
