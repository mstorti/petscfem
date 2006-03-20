// $Id: Mesh.cpp,v 1.1.2.4 2006/03/20 16:06:00 rodrigop Exp $

#include "Nodedata.h"
#include "Elemset.h"
#include "Mesh.h"

#include <fem.h>


PYPF_NAMESPACE_BEGIN


OptionTable*
Mesh::get_opt_table() const
{ 
  Mesh::Base* mesh = *this;
  if (mesh->global_options == NULL) mesh->global_options = new OptionTable;
  return mesh->global_options;
}


Mesh::~Mesh() 
{ 
  Mesh::Base* mesh = *this;
  /* nodedata */
  this->nodedata->decref();
  mesh->nodedata = NULL;
  /* elemset list */
  for (int i=0; i<this->elemsetlist.size(); i++)
    this->elemsetlist[i]->decref();
  PYPF_DELETE(da_destroy, mesh->elemsetlist);
  /* options */
  PYPF_DELETE_SCLR(mesh->global_options);
  /* base object */
  PYPF_DELETE_SCLR(mesh);
}

Mesh::Mesh() 
  : Ptr(new Mesh::Base), Object(),
    nodedata(new Nodedata), elemsetlist(0)
{ 
  Mesh::Base* mesh = *this;
  /* nodedata */
  this->nodedata->incref();
  mesh->nodedata = *(this->nodedata);
  /* elemset list */
  mesh->elemsetlist = da_create(sizeof(Elemset::Base*));
  /* options */
  mesh->global_options = new OptionTable;
}

Mesh::Mesh(const Mesh& _mesh) 
  : Ptr(new Mesh::Base), Object(_mesh),
    nodedata(_mesh.nodedata), elemsetlist(_mesh.elemsetlist)
{ 
  Mesh::Base* mesh = *this;
  /* nodedata */
  this->nodedata->incref();
  mesh->nodedata = *(this->nodedata);
  /* elemset list */
  mesh->elemsetlist = da_create_len(sizeof(Elemset::Base*),
				    this->elemsetlist.size());
  for (int i=0; i<this->elemsetlist.size(); i++) {
    Elemset* elemset = this->elemsetlist[i];
    elemset->incref();
    Elemset::Base* e = *elemset;
    da_set(mesh->elemsetlist, i, &e);
  }
  /* options */
  mesh->global_options = new OptionTable;
  this->setOptions(_mesh.getOptions());
}

Mesh::Mesh(Mesh::Base* _mesh)
  : Ptr(_mesh), Object(),
    nodedata(new Nodedata(_mesh->nodedata)), elemsetlist(0)
{  
  Mesh::Base* mesh = *this;
  /* nodedata */
  this->nodedata->incref();
  /* elemset list */
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
      elemset->incref();
      this->elemsetlist[i] = elemset;
    }
  }
  /* options */
  if (mesh->global_options == NULL) {
    mesh->global_options = new OptionTable;
  }
}

Mesh::Mesh(Nodedata* _nodedata,
	   const std::vector<Elemset*>& _elemsetlist) 
  : Ptr(new Mesh::Base), Object(),
    nodedata(_nodedata), elemsetlist(_elemsetlist)
{ 
  Mesh::Base* mesh = *this;
  /* nodedata */
  this->nodedata->incref();
  mesh->nodedata = *(this->nodedata);
  /* elemset list */
  mesh->elemsetlist = da_create_len(sizeof(Elemset::Base*),
				    this->elemsetlist.size());
  for (int i=0; i<this->elemsetlist.size(); i++) {
    Elemset* elemset = this->elemsetlist[i];
    elemset->incref();
    Elemset::Base* e = *elemset;
    da_set(mesh->elemsetlist, i, &e);
  }
  /* options */
  mesh->global_options = new OptionTable;
}


Nodedata*
Mesh::getNodedata() const
{
  return this->nodedata;
}

void
Mesh::setNodedata(Nodedata* nodedata)
{
  Mesh::Base* mesh = *this;
  nodedata->incref();
  this->nodedata->decref();
  this->nodedata = nodedata;
  mesh->nodedata = *nodedata;
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
  if (n == 0)      throw Error("empty elemset list");
  if (i<0 || i>=n) throw Error("index out of range");
  return this->elemsetlist[i];
}

void
Mesh::setElemset(int i, Elemset* elemset)
{
  int n = this->elemsetlist.size();
  if (n == 0)      throw Error("empty elemset list");
  if (i<0 || i>=n) throw Error("index out of range");
  /**/
  Mesh::Base* mesh = *this;
  elemset->incref();
  this->elemsetlist[i]->decref();
  this->elemsetlist[i] = elemset;
  Elemset::Base* e = *elemset;
  da_set(mesh->elemsetlist, i, &e);
}


void
Mesh::addElemset(Elemset* elemset)
{
  Mesh::Base* mesh = *this;
  elemset->incref();
  this->elemsetlist.push_back(elemset);
  Elemset::Base* e = *elemset;
  da_append(mesh->elemsetlist, &e);
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
  Mesh::Base* mesh = *this;
  /* nodedata */
  this->nodedata->decref();
  this->nodedata = new Nodedata;
  this->nodedata->incref();
  mesh->nodedata = *(this->nodedata);
  /* elemset list */
  for (int i=0; i<this->elemsetlist.size(); i++)
    this->elemsetlist[i]->decref();
  this->elemsetlist.resize(0);
  PYPF_DELETE(da_destroy, mesh->elemsetlist);
  mesh->elemsetlist = da_create(sizeof(Elemset::Base*));
  /* options */
  PYPF_DELETE_SCLR(mesh->global_options);
  mesh->global_options = new OptionTable;
}

void
Mesh::view() const
{
  printf("Nodedata\n");
  printf("--------\n");
  Nodedata::Base* nodedata = *(this->nodedata);
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
