// -*- mode: C++ -*- 
#ifndef NODEDATA_H
#define NODEDATA_H

#include <vector>
#include <string>

#include "fstack.h"
#include "texthash.h"
#include "getarray.h"

typedef string Field;
typedef int FieldIndx;
typedef list<string> SList;
typedef SMap map<string,int>;
typedef FieldM SMap;
typedef Slist FieldL;

class NodeData;

class FieldIndxV {
  vector<FieldIndx> v;
  int n;
public:
  friend NodeData;
}
  
void int2s(int n,string &s);

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** Contains (constant) data relative to nodes (may be coordinates and
    other). 
    @author M. Storti
    @param nnod number of nodes
    @param nu number of real quantities per node
*/ 
class NodeData {
  FieldL field_list,eq_list;
  FieldM fld_map,eq_map;
  double *data;
  TextHashTable *thash;
  int nnod,ndof,neqs;
  int part;
public:
  NodeData(FileStack &fstack);
  FieldIndx get_field_indx(const Field &field);
  get_field_indx(const FieldL &fieldv,FieldIndxV &field_indx_v);
  double val(const Field &field,const int node);
  void vals(const FieldIndxV &field_indx_v,int node,double *data);
  int fields() {return ndof;};
  int equations() {return neqs;};
  int size() {return nnod;};
#if 0
  void set_part_mode(int part_) {part=part_;};
  void set_proc(int node,int proc);
  int proc(int node);
#endif
};

#endif
