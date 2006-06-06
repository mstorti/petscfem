// $Id: Mesh.cpp,v 1.1.2.9 2006/06/06 15:44:27 dalcinl Exp $

#include "Mesh.h"

#include <fem.h>


PYPF_NAMESPACE_BEGIN


Mesh::~Mesh() 
{ 
  PYPF_DECREF(this->nodeset);
  for (int i=0; i<this->elemsetlist.size(); i++) {
    Elemset* elemset = this->elemsetlist[i];
    PYPF_DECREF(elemset);
  }
  /* base pointer */ Mesh::Base* mesh = *this;
  /*              */ if (mesh == NULL) return;
  /* options      */ mesh->global_options = NULL;
  /* nodedata     */ mesh->nodedata = NULL;
  /* elemset list */ PYPF_DELETE(da_destroy, mesh->elemsetlist);
  /* base object  */ delete mesh;
}

Mesh::Mesh(const Mesh& msh) 
  : Handle(new Mesh::Base), 
    Object(msh),
    nodeset(msh.nodeset), 
    elemsetlist(msh.elemsetlist)
{ 
  PYPF_INCREF(this->nodeset);
  for (int i=0; i<this->elemsetlist.size();
       this->elemsetlist[i++]->incref());
  /* base pointer */ Mesh::Base* mesh = *this;
  /* options      */ mesh->global_options = this->options;
  /* nodedata     */ mesh->nodedata = *(this->nodeset);
  /* elemset list */ mesh->elemsetlist = 
  /*              */    da_create_len(sizeof(Elemset::Base*),
  /*              */ 		      this->elemsetlist.size());
  /*              */ for (int i=0; i<this->elemsetlist.size(); i++) {
  /*              */   Elemset::Base* e = *this->elemsetlist[i];
  /*              */   da_set(mesh->elemsetlist, i, &e);
  /*              */ }
}

Mesh::Mesh(Nodeset& nodeset, const std::vector<Elemset*>& elemsetlist)
  : Handle(new Mesh::Base), 
    Object(nodeset.getComm()),
    nodeset(&nodeset),
    elemsetlist(elemsetlist)
{ 
  PYPF_INCREF(this->nodeset);
  for (int i=0; i<this->elemsetlist.size();
       this->elemsetlist[i++]->incref());
  /* base pointer */ Mesh::Base* mesh = *this;
  /* nodedata     */ mesh->nodedata = *(this->nodeset);
  /* elemset list */ mesh->elemsetlist = 
  /*              */   da_create_len(sizeof(Elemset::Base*),
  /*              */	             this->elemsetlist.size());
  /*              */ for (int i=0; i<this->elemsetlist.size(); i++) {
  /*              */    Elemset::Base* e = *this->elemsetlist[i];
  /*              */    da_set(mesh->elemsetlist, i, &e);
  /*              */ }
  /* options      */ mesh->global_options = this->options;
}

Nodeset&
Mesh::getNodeset() const
{
  return *this->nodeset;
}

// void
// Mesh::setNodeset(Nodeset& nodeset)
// {
//   PYPF_INCREF(&nodeset);
//   PYPF_DECREF(this->nodeset);
//   this->nodeset = &nodeset;
//   /* base pointer */ Mesh::Base* mesh = *this;
//   /* nodedata     */ mesh->nodedata = *(this->nodeset);
// }


int
Mesh::getSize() const
{
  return this->elemsetlist.size();
}
  
Elemset&
Mesh::getElemset(int i) const
{
  int n = this->elemsetlist.size();
  PYPF_ASSERT(i>=0 && i<n, "index out of range");
  return *this->elemsetlist[i];
}

// void
// Mesh::setElemset(int i, Elemset& elemset)
// {
//   int n = this->elemsetlist.size();
//   PYPF_ASSERT(i>=0 && i<n, "index out of range");

//   PYPF_INCREF(&elemset);
//   PYPF_DECREF(this->elemsetlist[i]);
//   this->elemsetlist[i] = &elemset;

//   /* base pointer */ Mesh::Base* mesh = *this;
//   /* base pointer */ Elemset::Base* e = elemset;
//   /* elemset list */ da_set(mesh->elemsetlist, i, &e);
// }

// void
// Mesh::delElemset(int i)
// {
//   int n = this->elemsetlist.size();
//   PYPF_ASSERT(i>=0 && i<n, "index out of range");
//   PYPF_ASSERT(0, "not implemented yet");
// }

// void
// Mesh::addElemset(Elemset& elemset)
// {
//   PYPF_INCREF(&elemset);
//   this->elemsetlist.push_back(&elemset);

//   /* base pointer */ Mesh::Base* mesh = *this;
//   /* base pointer */ Elemset::Base* e = elemset;
//   /* elemset list */ da_append(mesh->elemsetlist, &e);
// }

void
Mesh::view() const
{

  printf("Mesh Object:\n");

  printf("  Nodeset:\n");
  Nodeset::Base* nodedata = *(this->nodeset);
  printf("    nnod: %d, ndim: %d, nval=%d\n", 
	 nodedata->nnod, nodedata->ndim, nodedata->nu - nodedata->ndim);
  printf("  Elemset List\n");
  for (int i=0; i<this->elemsetlist.size(); i++) {
    Elemset::Base* elemset = *(this->elemsetlist[i]);
    printf("    %d -> type: %s, name: %s, sizes: (%d, %d)\n", 
	   i, elemset->type, elemset->name(), elemset->nelem, elemset->nel);
  }
}


PYPF_NAMESPACE_END
