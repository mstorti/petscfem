// -*- mode: C++ -*- 
#ifndef NODEDATA_H
#define NODEDATA_H

#include <vector>
#include <string>

typedef string Field;

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** Contains (constant) data relative to nodes (may be coordinates and
    other). 
    @author M. Storti
    @param nnod number of nodes
    @param nu number of real quantities per node
*/ 
class NodeData {
  vector<string> field_list;
public:
  int field_indx(Field &field)=0;
  double val(Field &field,int node)=0;
};

#endif
