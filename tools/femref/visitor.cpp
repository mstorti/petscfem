//__INSERT_LICENSE__
// $Id: visitor.cpp,v 1.10 2004/12/19 20:46:21 mstorti Exp $

#include <string>
#include <list>
#include "./hasher.h"

using namespace std;

#include "./femref.h"
#include "./gtemplates.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
UniformMesh::visitor::visitor() 
  : at_end(true), mesh(NULL), 
    etree_p(NULL), trace(0) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void UniformMesh::visitor::init(UniformMesh &mesh_a,int elem_a) {
  elem = elem_a;
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
  at_end = false;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool UniformMesh::visitor::so_next() {
  list<RefPathNode>::iterator ws, 
    w = ref_stack.begin();
  // here visit w...
  ElemRef::iterator q, qs, qrsib, qfather;
  q = w->splitter;
  int j = w->so_indx;
  if (trace) {
    printf("elem %d, (siblings ",elem);
    list<RefPathNode>::iterator 
      r = ref_stack.begin(), t; 
    t = r; t++;
    while (t!=ref_stack.end()) {
      printf("%d ",r->so_indx);
      r = t; t++;
    }
    printf("), ");
    w->go.print("go");
  }

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
    ws->splitter = qs;
    ws->go.make_canonical();
    ws->so_indx = 0;
    return true;
  } else {
    // `q' is a leave for GO's (sure it isn't a
    // regular node for splitters). Try to find a
    // right sibling, or a father that has a right
    // sibling
    while (true) {
      // Check if we are at the root
      if (ref_stack.size()<=1) {
	at_end = true;
	return false;
      }
      
      ws = ref_stack.begin(); 
      j = ws->so_indx;
      w=ws; w++;
      qfather = w->splitter;
      assert(qfather != etree_p->end());
      const Splitter *s = qfather->splitter;

      ref_stack.pop_front();

      if (j<s->size()-1) {
	int jsib = j+1;
	qrsib = ws->splitter;
	ref_stack.push_front(RefPathNode());
	ws = ref_stack.begin();
	w = ws; w++;
	// Build `ws' from GeomObject `w' (parent) and splitter `s'
	// and subobject index `jsib'
	mesh->set(w->go,s,jsib,ws->go);
	ws->go.make_canonical();

	// Find next node on the splitting tree or end()
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
init(UniformMesh &mesh_a) { 
  init(mesh_a,0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool UniformMesh::visitor::
next() { 
  if (so_next()) return true;
  elem++;
  if (elem >= mesh->nelem) return false;
  init(*mesh,elem);
  return true;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool UniformMesh::visitor::
so_end() { return at_end; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool UniformMesh::visitor::
end() { return at_end && elem>= mesh->nelem; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void UniformMesh::visitor::
refine(const Splitter* s) {
  assert(is_leave());
  list<RefPathNode>::iterator 
    w = ref_stack.begin();
  ElemRef::iterator 
    &q = w->splitter;
  q = etree_p->insert(q,ElemRefNode());
  int j = w->so_indx;
  // The splitter should be returned
  // by the refinement function
  q->splitter = s;
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool UniformMesh::visitor::level_so_next() {
  list<RefPathNode>::iterator ws, 
    w = ref_stack.begin();
  // here visit w...
  ElemRef::iterator q, qs, qrsib, qfather;
  q = w->splitter;
  int j = w->so_indx;
  if (trace) {
    printf("elem %d, (siblings ",elem);
    list<RefPathNode>::iterator 
      r = ref_stack.begin(), t; 
    t = r; t++;
    while (t!=ref_stack.end()) {
      printf("%d ",r->so_indx);
      r = t; t++;
    }
    printf("), ");
    w->go.print("go");
  }

  // Doesn't try to follow the childs. Go
  // same level or higher. Try to find a
  // right sibling, or a father that has a right
  // sibling
  while (true) {
    // Check if we are at the root
    if (ref_stack.size()<=1) {
      at_end = true;
      return false;
    }
    
    ws = ref_stack.begin(); 
    j = ws->so_indx;
    w=ws; w++;
    qfather = w->splitter;
    assert(qfather != etree_p->end());
    const Splitter *s = qfather->splitter;
    
    ref_stack.pop_front();
    
    if (j<s->size()-1) {
      int jsib = j+1;
      qrsib = ws->splitter;
      ref_stack.push_front(RefPathNode());
      ws = ref_stack.begin();
      w = ws; w++;
      // Build `ws' from GeomObject `w' (parent) and splitter `s'
      // and subobject index `jsib'
      mesh->set(w->go,s,jsib,ws->go);
      ws->go.make_canonical();
      
      // Find next node on the splitting tree or end()
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool UniformMesh::visitor::
level_next() { 
  if (level_so_next()) return true;
  elem++;
  if (elem >= mesh->nelem) return false;
  init(*mesh,elem);
  return true;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool UniformMesh::visitor::
next(int level) {
  if (ref_level()<level) return next();
  else return level_next();
}
