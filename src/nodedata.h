// -*- mode: C++ -*- 
#ifndef NODEDATA_H
#define NODEDATA_H

class Field {
public:
  virtual boolean operator=(Field &)=0;
}

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** Contains (constant) data relative to nodes (may be coordinates and
    other). 
    @author M. Storti
    @param nnod number of nodes
    @param nu number of real quantities per node
*/ 
class NodeData {
public:
  virtual int field_indx(Field &field)=0;
  virtual double val(Field &field,int node)=0;
};



#endif
