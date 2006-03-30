// $Id: Mesh.cpp,v 1.1.2.6 2006/03/30 15:40:05 rodrigop Exp $

#include "Nodeset.h"
#include "Elemset.h"
#include "Mesh.h"

#include <fem.h>


PYPF_NAMESPACE_BEGIN


Mesh::~Mesh() 
{ 
  this->nodedata->decref();
  for (int i=0; i<this->elemsetlist.size();
       this->elemsetlist[i++]->decref());
  /* base pointer */ Mesh::Base* mesh = *this;
  /* options      */ mesh->global_options = NULL;
  /* nodedata     */ mesh->nodedata = NULL;
  /* elemset list */ PYPF_DELETE(da_destroy, mesh->elemsetlist);
  /* base object  */ delete mesh;
}

Mesh::Mesh() 
  : Handle(new Mesh::Base), Object(),
    nodedata(new Nodeset), elemsetlist(0)
{ 
  this->nodedata->incref();
  /* base pointer */ Mesh::Base* mesh = *this;
  /* options      */ mesh->global_options = this->options;
  /* nodedata     */ mesh->nodedata = *(this->nodedata); 
  /* elemset list */ mesh->elemsetlist = da_create(sizeof(Elemset::Base*));
}

Mesh::Mesh(const Mesh& msh) 
  : Handle(new Mesh::Base), Object(msh),
    nodedata(msh.nodedata), elemsetlist(msh.elemsetlist)
{ 
  this->nodedata->incref();
  for (int i=0; i<this->elemsetlist.size(); 
       this->elemsetlist[i++]->incref());
  /* base pointer */ Mesh::Base* mesh = *this;
  /* options      */ mesh->global_options = this->options;
  /* nodedata     */ mesh->nodedata = *(this->nodedata);
  /* elemset list */ mesh->elemsetlist = 
  /*              */    da_create_len(sizeof(Elemset::Base*),
  /*              */ 		      this->elemsetlist.size());
  /*              */ for (int i=0; i<this->elemsetlist.size(); i++) {
  /*              */   Elemset::Base* e = *this->elemsetlist[i];
  /*              */   da_set(mesh->elemsetlist, i, &e);
  /*              */ }
}

Mesh::Mesh(Mesh::Base* msh)
  : Handle(msh), Object(),
    nodedata(new Nodeset(msh->nodedata)), elemsetlist(0)
{  
  Mesh::Base* mesh = *this;
  if (mesh->elemsetlist == NULL) {
    mesh->elemsetlist = da_create(sizeof(Elemset::Base*));
  }
  else {
    int size = int(da_length(mesh->elemsetlist));
    this->elemsetlist.resize(size);
    for (int i=0; i<size; i++) {
      Elemset::Base** e;
      e = reinterpret_cast<Elemset::Base**>(da_ref(mesh->elemsetlist, i));
      Elemset* elemset = new Elemset(*e);
      this->elemsetlist[i] = elemset;
    }
  }
  if (mesh->global_options == NULL) 
    mesh->global_options = this->options;
  else
    this->options = mesh->global_options;

  this->nodedata->incref();
  for (int i=0; i<this->elemsetlist.size(); 
       this->elemsetlist[i++]->incref());
}

Mesh::Mesh(Nodeset* _nodedata,
	   const std::vector<Elemset*>& _elemsetlist) 
  : Handle(new Mesh::Base), Object(),
    nodedata(_nodedata), elemsetlist(_elemsetlist)
{ 
  this->nodedata->incref();
  for (int i=0; i<this->elemsetlist.size(); 
       this->elemsetlist[i++]->incref());
  /* base pointer */ Mesh::Base* mesh = *this;
  /* nodedata     */ mesh->nodedata = *(this->nodedata);
  /* elemset list */ mesh->elemsetlist = 
  /*              */   da_create_len(sizeof(Elemset::Base*),
  /*              */	             this->elemsetlist.size());
  /*              */ for (int i=0; i<this->elemsetlist.size(); i++) {
  /*              */    Elemset::Base* e = *this->elemsetlist[i];
  /*              */    da_set(mesh->elemsetlist, i, &e);
  /*              */ }
  /* options      */ mesh->global_options = this->options;
}


Nodeset*
Mesh::getNodeset() const
{
  return this->nodedata;
}

void
Mesh::setNodeset(Nodeset* nodedata)
{
  nodedata->incref();
  this->nodedata->decref();
  this->nodedata = nodedata;
  /* base pointer */ Mesh::Base* mesh = *this;
  /* nodedata     */ mesh->nodedata = *(this->nodedata);
}


int
Mesh::getSize() const
{
  return this->elemsetlist.size();
}
  
Elemset*
Mesh::getElemset(int i) const
{
  int n = this->elemsetlist.size();
  if (n==0)      throw Error("empty elemset list");
  if (i<0||i>=n) throw Error("index out of range");
  return this->elemsetlist[i];
}

void
Mesh::setElemset(int i, Elemset* elemset)
{
  int n = this->elemsetlist.size();
  if (n==0)      throw Error("empty elemset list");
  if (i<0||i>=n) throw Error("index out of range");

  elemset->incref();
  this->elemsetlist[i]->decref();
  this->elemsetlist[i] = elemset;

  /* base pointer */ Mesh::Base* mesh = *this;
  /* base pointer */ Elemset::Base* e = *elemset;
  /* elemset list */ da_set(mesh->elemsetlist, i, &e);
}


void
Mesh::addElemset(Elemset* elemset)
{
  elemset->incref();
  this->elemsetlist.push_back(elemset);

  /* base pointer */ Mesh::Base* mesh = *this;
  /* base pointer */ Elemset::Base* e = *elemset;
  /* elemset list */ da_append(mesh->elemsetlist, &e);
}

void
Mesh::setUp()
{
  this->nodedata->setUp();
  for (int i=0; i<this->elemsetlist.size(); i++)
    this->elemsetlist[i]->setUp();
}

void
Mesh::clear()
{
  this->options.clear();
  this->nodedata->decref();
  this->nodedata = new Nodeset;
  this->nodedata->incref();
  for (int i=0; i<this->elemsetlist.size(); 
       this->elemsetlist[i++]->decref());
  this->elemsetlist.resize(0);

  /* base pointer */ Mesh::Base* mesh = *this;
  /* options      */ mesh->global_options = this->options;
  /* nodedata     */ mesh->nodedata = *(this->nodedata);
  /* elemset list */ PYPF_DELETE(da_destroy, mesh->elemsetlist);
  /*              */ mesh->elemsetlist = da_create(sizeof(Elemset::Base*));

}

void
Mesh::view() const
{
  printf("Nodeset\n");
  printf("--------\n");
  Nodeset::Base* nodedata = *(this->nodedata);
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
