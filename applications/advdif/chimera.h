#ifndef PETSCFEM_CHIMERA_H
#define PETSCFEM_CHIMERA_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Marks the nodes that are bdry (internal and external)
class chimera_hook_t : public Hook {
public:
  virtual void mark_bdry_nodes(set<int> &ebdry,set<int> &ibdry) { }
};

#endif
