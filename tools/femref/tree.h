// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: tree.h,v 1.4 2004/12/26 03:09:24 mstorti Exp $
#ifndef PETSCFEM_TREE_H
#define PETSCFEM_TREE_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
class tree {
public:
  class iterator;
  class cell {
    friend class tree;
    friend class iterator;
    T t;
    cell *right, *left_child, *father;
    cell() : right(NULL), left_child(NULL) {}
  public:
    tree<T> *tree_owner() {
      cell *q = this;
      while (q->father) q = q->father;
      return (tree<T>*)q->right;
    }
  };
private:
  cell *header;

  iterator tree_copy_aux(iterator nq,
			 tree<T> &TT,iterator nt) {
    nq = insert(nq,*nt);
    iterator
      ct = nt.lchild(),
      cq = nq.lchild();
    while (ct!=TT.end()) {
      cq = tree_copy_aux(cq,TT,ct);
      ct = ct.right();
      cq = cq.right();
    }
  }
 public:
  static int cell_count_m;
  static int cell_count() { return cell_count_m; }
  class iterator {
  private:
    friend class tree;
    tree<T> *tree_p;
    cell *ptr,*prev,*father_p;
    iterator(cell *p,cell *prev_a,cell *father_a,tree<T> *t) : ptr(p), 
      prev(prev_a), tree_p(t), father_p(father_a) { }
  public:
    cell *cell_ptr() { return ptr; }
    iterator(const iterator &q) {
      ptr = q.ptr;
      prev = q.prev; 
      father_p =q.father_p;
      tree_p = q.tree_p;
    }
    T &operator*() { return ptr->t; }
    T *operator->() { return &ptr->t; }
    bool operator!=(iterator q) { return ptr!=q.ptr; }
    bool operator==(iterator q) { return ptr==q.ptr; }
    iterator() : ptr(NULL), prev(NULL), 
		 tree_p(NULL), father_p(NULL) { }

    iterator lchild() { 
      return iterator(ptr->left_child,
		      NULL,ptr,tree_p); 
    }
    iterator right() { 
      return iterator(ptr->right,ptr,
		      father_p,tree_p); 
    }
    iterator father() { 
      cell *father_p = ptr->father;
      if (father_p == tree_p->header) return tree_p->end();
      cell *grand_father = father_p->father;
      cell *father_prev = grand_father->left_child;
      assert(father_prev);
      if(father_prev == father_p) 
	return iterator(father_p,NULL,grand_father,tree_p);
      while (father_prev 
	     && father_prev->right!=father_p)
	father_prev = father_prev->right;
      return iterator(father_p,father_prev,grand_father,tree_p);
    }

    // Prefix:
    iterator operator++() {
      *this = right();
      return *this;
    }
    // Postfix:
    iterator operator++(int) {
      iterator q = *this;
      *this = right();
      return q;
    }
    tree<T> *tree_owner() {
      return tree_p;
    }
  };

  tree() {
    header = new cell;
    cell_count_m++;
    // This is somewhat tricky. We store in the right
    // field of the header cell the pointer to the
    // tree. So that from any cell we can retrive
    // the whole tree.
    // The tree pointer is converted to a cell pointer,
    // but I think this is OK. 
    header->right = (cell *)this;
    header->left_child = NULL;
    header->father = NULL;
  }
  tree<T>(const tree<T> &TT) { 
    if (&TT != this) {
      header = new cell;
      cell_count_m++;
      header->right = NULL;
      header->left_child = NULL;
      header->father = NULL;
      tree<T> &TTT = (tree<T> &) TT;
      if (TTT.begin()!=TTT.end()) 
	tree_copy_aux(begin(),TTT,TTT.begin()); 
    }
  }
  ~tree() { clear(); delete header; cell_count_m--; }
  iterator insert(iterator p,T t) {
    assert(!(p.ptr && p.ptr->father == header));
    cell *c = new cell;
    cell_count_m++;
    c->right = p.ptr;
    c->t = t;
    c->father = p.father_p;
    p.ptr = c;
    if (p.prev) p.prev->right = c;
    else p.father_p->left_child = c;	
    return p;
  }
  iterator erase(iterator p) {
    if(p==end()) return p;
    iterator c = p.lchild();
    while (c!=end()) c = erase(c);
    cell *q = p.ptr;
    p.ptr = p.ptr->right;
    if (p.prev) p.prev->right = p.ptr;
    else p.father_p->left_child = p.ptr;
    delete q;
    cell_count_m--;
    return p;
  }

  iterator splice(iterator to,iterator from) {
    assert(!(to.father_p==header && to.ptr));
    cell *c = from.ptr;

    if (from.prev) from.prev->right = c->right;
    else from.father_p->left_child = c->right;
    c->father = to.father_p;

    c->right = to.ptr;
    to.ptr = c;
    if (to.prev) to.prev->right = c;
    else to.father_p->left_child = c;	

    return to;
  }
  iterator find(T t) { return find(t,begin()); }
  iterator find(T t,iterator p) {
    if(p==end() || p.ptr->t == t) return p;
    iterator q,c = p.lchild();
    while (c!=end()) {
      q = find(t,c);
      if (q!=end()) return q;
      else c++;
    }
    return iterator();
  }
  void clear() { erase(begin()); }
  iterator begin() { 
    return iterator(header->left_child,NULL,header,this); 
  }
  iterator end() { return iterator(); }

  //PP>if 0
  void print_prev(iterator p) { 
    if (p==end()) return;
    cout << "(" << p.ptr << "," << p.ptr->t << ")" << endl;
    iterator c = p.lchild();
    while (c!=end()) print_prev(c++);
  }
  void print_prev() { print_prev(begin()); }

  void print_post(iterator p) { 
    if (p==end()) return;
    iterator c = p.lchild();
    while (c!=end()) print_post(c++);
    cout << "(" << p.ptr << "," << p.ptr->t << ")" << endl;
  }
  void print_post() { print_post(begin()); }

  void lisp_print(iterator n) {
    if (n==end()) return;
    iterator c = n.lchild();
    bool is_leaf = c==end();
    if (is_leaf) cout << *n;
    else {
      cout << "(" << *n;
      while (c!=end()) {
	cout << " ";
	lisp_print(c++);
      }
      cout << ")";
    }
  }
  void lisp_print() { lisp_print(begin()); }
};

template<class T>
int tree<T>::cell_count_m = 0;

#endif
