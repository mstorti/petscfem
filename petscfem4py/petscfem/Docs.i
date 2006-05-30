// -*- c++ -*-
// $Id: Docs.i,v 1.1.2.4 2006/05/30 20:20:25 dalcinl Exp $

%define %docstring(node,...) 
#if #__VA_ARGS__ == ""
  %feature("autodoc", 2) node;
#else
  %feature("autodoc", __VA_ARGS__) node;
#endif
%enddef

%define %doctypemap(Type)
%typemap(doc) Type, Type*, Type&, const Type&
  "$1_name: Type instance";
%enddef

namespace std { }
%typemap(doc, type="string") 
  std::string, const std::string&
  "$1_name: string value";

%typemap(doc, type="map<string,string>") 
  std::map<std::string,std::string>,
  const std::map<std::string,std::string>&
  "$1_name: dict instance";


%import namespace.h


PYPF_NAMESPACE_BEGIN

%feature("autodoc");

%doctypemap(Comm);
%docstring(Comm);
%docstring(Comm::Comm);
%docstring(Comm::getComm);
%docstring(Comm::setComm);
%docstring(Comm::operator==);
%docstring(Comm::operator!=);
%docstring(Comm::getSize);
%docstring(Comm::getRank);


%doctypemap(Object);
%docstring(Object);
%docstring(Object::getComm);
%docstring(Object::setComm);
%docstring(Object::hasOption);
%docstring(Object::getOption);
%docstring(Object::setOption);
%docstring(Object::getOptions, "getOptions(self) -> dict");
%docstring(Object::setOptions);
%docstring(Object::addOptions);
%docstring(Object::delOptions);
%docstring(Object::operator==);
%docstring(Object::operator!=);

%typemap(doc,name="node",type="double array[]") 
  (int n, const double data[])
  "node: double array[] value";
%typemap(doc, name="nodedata",type="double array[]") 
  (int nnod, int nval, const double data[])
  "nodedata: double array[] value";

%doctypemap(Nodeset);
%docstring(Nodeset);
%docstring(Nodeset::Nodeset);
%docstring(Nodeset::getData,     "getData(self, int i) -> double array[]") ;
%docstring(Nodeset::setData);
%docstring(Nodeset::getDataSize, "getDataSize(self) -> (int, int)") ;
%docstring(Nodeset::getNode,     "getNode(self, int i) -> double array[]");
%docstring(Nodeset::setNode);
%docstring(Nodeset::getSize);
%docstring(Nodeset::getDim);
%docstring(Nodeset::setDim);
%docstring(Nodeset::clear);
%docstring(Nodeset::view);
%docstring(Nodeset::__len__);
%docstring(Nodeset::__getitem__, "__getitem__(self, int i) -> double array[]");
%docstring(Nodeset::__setitem__);
%docstring(Nodeset::__iter__,    "__iter__(self) -> iterator");


%typemap(doc,name="elem",type="int array[]")
  (int n, const int elem[])
  "elem: int array[] value";
%typemap(doc,name="elemdata",type="int array[]")
  (int nelem, int nel, const int icone[])
  "elemdata: double array[] value";

%doctypemap(Elemset);
%docstring(Elemset);
%docstring(Elemset::Elemset);
%docstring(Elemset::getType);
%docstring(Elemset::setType);
%docstring(Elemset::getName);
%docstring(Elemset::setName);
%docstring(Elemset::getData,     "getData(self) -> int array[]");
%docstring(Elemset::setData);
%docstring(Elemset::getDataSize, "getDataSize(self) -> (int, int)");
%docstring(Elemset::getElem,     "getElem(self, int i) -> int array[]") ;
%docstring(Elemset::setElem);
%docstring(Elemset::getSize);
%docstring(Elemset::clear);
%docstring(Elemset::view);
%docstring(Elemset::__len__);
%docstring(Elemset::__getitem__, "__getitem__(self, int i) -> int array[]");
%docstring(Elemset::__setitem__);
%docstring(Elemset::__iter__,    "__iter__(self) -> iterator");

%typemap(doc, type="list<Elemset>")
  std::vector<PYPF_NAMESPACE::Elemset*>&,
  const std::vector<PYPF_NAMESPACE::Elemset*>&
  "$1_name: list of Elemset instances";

%doctypemap(Mesh);
%docstring(Mesh);
%docstring(Mesh::Mesh);
%docstring(Mesh::getNodeset);
%docstring(Mesh::setNodeset);
%docstring(Mesh::getElemset);
%docstring(Mesh::setElemset);
%docstring(Mesh::delElemset);
%docstring(Mesh::addElemset);
%docstring(Mesh::getSize);
%docstring(Mesh::clear);
%docstring(Mesh::view);
%docstring(Mesh::__len__);
%docstring(Mesh::__getitem__);
%docstring(Mesh::__setitem__);
%docstring(Mesh::__delitem__);
%docstring(Mesh::__iter__, "__iter__(self) -> iterator");


%doctypemap(Amplitude);
%docstring(Amplitude);
%docstring(Amplitude::Amplitude);
%docstring(Amplitude::operator());


%typemap(doc,name="node, field, value",type="int[], int[], double[]")
  (int n, const int node[],const int field[],const double value[])
  "node, field, value: int[], int[], double[] array values";
%typemap(doc,name="node, field, coeff",type="int[], int[], double[]")
  (int n, const int node[],const int field[],const double coeff[])
  "node, field, coeff: int[], int[], double[] array values";

%doctypemap(Dofset);
%docstring(Dofset);
%docstring(Dofset::Dofset);
%docstring(Dofset::addFixations);
%docstring(Dofset::addConstraints);
%docstring(Dofset::clear);


%doctypemap(DofMap);
%docstring(DofMap);
%docstring(DofMap::DofMap);
%docstring(DofMap::getSize);
%docstring(DofMap::getLocalSize);
%docstring(DofMap::getSizes, "getSizes(self) -> (int, int)");
%docstring(DofMap::getRange, "getRange(self) -> (int, int)");
%docstring(DofMap::getDist,  "getDist(self) -> int array[]");
%docstring(DofMap::getNNod);
%docstring(DofMap::getNDof);
%docstring(DofMap::apply);
%docstring(DofMap::solve);
%docstring(DofMap::view);


%doctypemap(Domain);
%docstring(Domain);
%docstring(Domain::Domain);
%docstring(Domain::getMesh);
%docstring(Domain::getDofMap);
%docstring(Domain::getDim);
%docstring(Domain::getSize);
%docstring(Domain::getSizes,    "getSizes(self) -> (int, int)");
%docstring(Domain::getDofSize);
%docstring(Domain::getDofSizes, "getDofSizes(self) -> (int, int)");
%docstring(Domain::getDofRange, "getDofRange(self) -> (int, int)");
%docstring(Domain::getDofDist,  "getDofDist(self) -> int array[]");

%docstring(Domain::getOwnedDofs, "getOwnedDofs(self) -> (int, int)");
%docstring(Domain::getGhostDofs, "getGhostDofs(self) -> int array[]");
%docstring(Domain::getLocalDofs, "getLocalDofs(self) -> int array[]");

%typemap(doc,name="state",type="double[]")
  (int ns, const double state[])
  "state: double[] array";
%typemap(doc,name="nodes",type="int[]")
  (int nn, const int nodes[])
  "nodes: int[] array";
%typemap(doc,name="fields",type="int[]")
  (int nf, const int fields[])
  "fields: int[] array";
%typemap(doc,name="values",type="double[]")
  std::vector<double>& values,
  (std::pair<int,int>& shape, std::vector<double>& values)
  "values: double[] array (output)";

%doctypemap(Application);
%docstring(Application);
%docstring(Application::Application);
%docstring(Application::getDomain);
%docstring(Application::getNodalValues);


%doctypemap(NvrStks);
%docstring(NvrStks);
%docstring(NvrStks::NvrStks);
%docstring(NvrStks::getAlpha);
%docstring(NvrStks::getSteady);
%docstring(NvrStks::setAlpha);
%docstring(NvrStks::setSteady);
%docstring(NvrStks::assemble);


PYPF_NAMESPACE_END
