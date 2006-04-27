// $Id: Mesh.cpp,v 1.1.2.7 2006/04/27 19:09:17 rodrigop Exp $

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
  /* options      */ mesh->global_options = NULL;
  /* nodedata     */ mesh->nodedata = NULL;
  /* elemset list */ PYPF_DELETE(da_destroy, mesh->elemsetlist);
  /* base object  */ delete mesh;
}

Mesh::Mesh() 
  : Handle(new Mesh::Base), 
    Object(),
    nodeset(new Nodeset), 
    elemsetlist()
{ 
  PYPF_INCREF(this->nodeset);
  /* base pointer */ Mesh::Base* mesh = *this;
  /* options      */ mesh->global_options = this->options;
  /* nodedata     */ mesh->nodedata = *(this->nodeset); 
  /* elemset list */ mesh->elemsetlist = da_create(sizeof(Elemset::Base*));
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

Mesh::Mesh(Nodeset& nodeset,
	   const std::vector<Elemset*>& elemsetlist)
  : Handle(new Mesh::Base), 
    Object(),
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

void
Mesh::setNodeset(Nodeset& nodeset)
{
  PYPF_INCREF(&nodeset);
  PYPF_DECREF(this->nodeset);
  this->nodeset = &nodeset;
  /* base pointer */ Mesh::Base* mesh = *this;
  /* nodedata     */ mesh->nodedata = *(this->nodeset);
}


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

void
Mesh::setElemset(int i, Elemset& elemset)
{
  int n = this->elemsetlist.size();
  PYPF_ASSERT(i>=0 && i<n, "index out of range");

  PYPF_INCREF(&elemset);
  PYPF_DECREF(this->elemsetlist[i]);
  this->elemsetlist[i] = &elemset;

  /* base pointer */ Mesh::Base* mesh = *this;
  /* base pointer */ Elemset::Base* e = elemset;
  /* elemset list */ da_set(mesh->elemsetlist, i, &e);
}

void
Mesh::delElemset(int i)
{
  int n = this->elemsetlist.size();
  PYPF_ASSERT(i>=0 && i<n, "index out of range");
  PYPF_ASSERT(0, "not implemented yet");
}

void
Mesh::addElemset(Elemset& elemset)
{
  PYPF_INCREF(&elemset);
  this->elemsetlist.push_back(&elemset);

  /* base pointer */ Mesh::Base* mesh = *this;
  /* base pointer */ Elemset::Base* e = elemset;
  /* elemset list */ da_append(mesh->elemsetlist, &e);
}

void
Mesh::clear()
{
  this->options.clear();
  PYPF_DECREF(this->nodeset);
  this->nodeset = new Nodeset;
  PYPF_INCREF(this->nodeset);
  for (int i=0; i<this->elemsetlist.size(); 
       this->elemsetlist[i++]->decref());
  this->elemsetlist.resize(0);

  /* base pointer */ Mesh::Base* mesh = *this;
  /* options      */ mesh->global_options = this->options;
  /* nodedata     */ mesh->nodedata = *(this->nodeset);
  /* elemset list */ PYPF_DELETE(da_destroy, mesh->elemsetlist);
  /*              */ mesh->elemsetlist = da_create(sizeof(Elemset::Base*));

}

void
Mesh::view() const
{
  printf("Nodeset\n");
  printf("--------\n");
  Nodeset::Base* nodedata = *(this->nodeset);
  printf("nnod: %d, ndim: %d\n", nodedata->nnod, nodedata->ndim);
  printf("\n");


  printf("Elemset List\n");
  printf("------------\n");
  for (int i=0; i<this->elemsetlist.size(); i++) {
    Elemset::Base* elemset = *(this->elemsetlist[i]);
    printf("Elemset %d: type: %s, name: %s, size: (%d, %d)\n", i, 
	   elemset->type, elemset->name(),
	   elemset->nelem, elemset->nel);
  }
}


PYPF_NAMESPACE_END
