//__INSERT_LICENSE__
// $Id: combiner.cpp,v 1.1 2004/12/25 03:45:04 mstorti Exp $

#include <string>
#include <list>
#include <limits.h>
#include "./hasher.h"

using namespace std;

#include "./femref.h"
#include "./gtemplates.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void NodeInfoRef::print() {
  printf("(coords \n");
  for (int j=0; j<coords.size(); j++) 
    printf("%f ",coords[j]);
  printf("\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void 
LinearCombiner::
combine(int tag,
	int n,const int *nodes,
	int new_node,
	NodeInfoMapT &node_info_map) {
  const int &ndim = mesh->ndim;
  if (n==0 ) {
    printf("node-comb (%d)\n",new_node);
    if (node_info_map.find(new_node) 
	== node_info_map.end()) {
      NodeInfoRef *ni_p = new NodeInfoRef;
      vector<double> &coords = ni_p->coords;
      coords.resize(ndim);
      for (int j=0; j<ndim; j++) 
	coords[j] = mesh->coords.e(new_node,j);
      node_info_map[new_node] = ni_p;
    } 
  } else if (n==2) {
    printf("node-comb (%d,%d) -> %d\n",
	   nodes[0],nodes[1],new_node);
    
    NodeInfoMapT::iterator q;
      
    q = node_info_map.find(nodes[0]);
    if(q == node_info_map.end()) {
      printf("slot for node %d in node_info_map is empty\n",
	     nodes[0]);
      assert(0);
    }
    const NodeInfoRef *ni1_p 
      = dynamic_cast<const NodeInfoRef *>(q->second);
      
    q = node_info_map.find(nodes[1]);
    if(q == node_info_map.end()) {
      printf("slot for node %d in node_info_map is empty\n",
	     nodes[1]);
      assert(0);
    }
    const NodeInfoRef *ni2_p 
      = dynamic_cast<const NodeInfoRef *>(q->second);
      
    q = node_info_map.find(new_node);
    if(q != node_info_map.end()) {
      printf("slot for node %d in node_info_map is already filled\n",
	     new_node);
      assert(0);
    }
    NodeInfoRef *ni_p = new NodeInfoRef;
    ni_p->coords.resize(ndim);
    for (int j=0; j<ndim; j++) {
      ni_p->coords[j] = 0.5*(ni1_p->coords[j]
			     +ni2_p->coords[j]);
    }
    node_info_map[new_node] = ni_p;
  } else assert(0);
}

LinearCombiner linear_combiner;
