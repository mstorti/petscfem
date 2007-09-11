// -*- c++ -*-
// $Id$

%define %doctypemap(Class)
%typemap(doc) Class, Class*, Class&, const Class&
  "$1_name: Class instance";
%enddef

%define %docstring(node,...) 
#if #__VA_ARGS__ == ""
  %feature("autodoc", 2) node;
#else
  %feature("autodoc", __VA_ARGS__) node;
#endif
%enddef

namespace std { }

%typemap(doc, type="string") 
  std::string, const std::string&
  "$1_name: string value";

%typemap(doc, type="pair<int,int>") 
  std::pair<int,int>, const std::pair<int,int>&
  "$1_name: pair<int,int> value";

%typemap(doc, type="map<string,string>") 
  std::map<std::string,std::string>, const std::map<std::string,std::string>&
  "$1_name: dict instance";


PF4PY_NAMESPACE_BEGIN
%feature("autodoc");
PF4PY_NAMESPACE_END


PF4PY_NAMESPACE_BEGIN

%doctypemap(Comm);
%docstring(Comm);
%docstring(Comm::Comm);
%docstring(Comm::operator==);
%docstring(Comm::operator!=);
%docstring(Comm::getSize);
%docstring(Comm::getRank);

%doctypemap(Object);
%docstring(Object);
%docstring(Object::Object);
%docstring(Object::operator==);
%docstring(Object::operator!=);

%doctypemap(DTable<int>);
%docstring(DTable<int>);
%docstring(DTable<int>::DTable);
%docstring(DTable<int>::getSize)
%docstring(DTable<int>::getShape, "getShape(self) -> (int, int)")
%docstring(DTable<int>::getArray);

%doctypemap(DTable<double>);
%docstring(DTable<double>);
%docstring(DTable<double>::DTable);
%docstring(DTable<double>::getSize)
%docstring(DTable<double>::getShape, "getShape(self) -> (int, int)")
%docstring(DTable<double>::getArray);

%doctypemap(Options);
%docstring(Options);
%docstring(Options::Options);
%docstring(Options::has);
%docstring(Options::get);
%docstring(Options::set);
%docstring(Options::del);
%docstring(Options::add);
%docstring(Options::update);
%docstring(Options::size);
%docstring(Options::clear);
%docstring(Options::todict, "todict(self) -> dict");

// --

%doctypemap(Elemset);
%docstring(Elemset);
%docstring(Elemset::Elemset);
%docstring(Elemset::getType);
%docstring(Elemset::setType);
%docstring(Elemset::getData);
%docstring(Elemset::setData);
%docstring(Elemset::getOptions);
%docstring(Elemset::setOptions);

%doctypemap(Mesh);
%docstring(Mesh);
%docstring(Mesh::Mesh);
%docstring(Mesh::getNodedata);
%docstring(Mesh::setNodedata);
%docstring(Mesh::getField);
%docstring(Mesh::setField);
%docstring(Mesh::getSize);
%docstring(Mesh::getElemset);
%docstring(Mesh::addElemset);
%docstring(Mesh::getOptions);
%docstring(Mesh::setOptions);

// --

%doctypemap(Amplitude);
%docstring(Amplitude);
%docstring(Amplitude::Amplitude);
%docstring(Amplitude::operator());

%doctypemap(Dofset);
%docstring(Dofset);
%docstring(Dofset::Dofset);

%docstring(Dofset::getComm);
%docstring(Dofset::getNNod);
%docstring(Dofset::getNDof);

%docstring(Dofset::getWeights);
%docstring(Dofset::setWeights);
%docstring(Dofset::getSizes);
%docstring(Dofset::getRange);
%docstring(Dofset::getDist);

%docstring(Dofset::getGhostDofs);
%docstring(Dofset::getLocalDofs);
%docstring(Dofset::getFieldDofs);

// --

%doctypemap(AppCtx);
%docstring(AppCtx);
%docstring(AppCtx::AppCtx);
%docstring(AppCtx::getAlpha);
%docstring(AppCtx::setAlpha);
%docstring(AppCtx::getSteady);
%docstring(AppCtx::setSteady);

// --

%doctypemap(Domain);
%docstring(Domain);
%docstring(Domain::Domain);

%docstring(Domain::getComm);
%docstring(Domain::getNDim);
%docstring(Domain::getNNod);
%docstring(Domain::getNDof);

%docstring(Domain::getType);

%docstring(Domain::getOptions);
%docstring(Domain::setOptions);

%docstring(Domain::getNodedata);
%docstring(Domain::setNodedata);

%docstring(Domain::getField);
%docstring(Domain::setField);

%docstring(Domain::getElemset);
%docstring(Domain::setElemset);
%docstring(Domain::addElemset);

%docstring(Domain::setFixation);
%docstring(Domain::setPeriodic);
%docstring(Domain::setConstraint);

%docstring(Domain::getMesh);
%docstring(Domain::getDofset);
%docstring(Domain::getAppCtx);

%docstring(Domain::setUp);

%docstring(Domain::allocateState);
%docstring(Domain::allocateResidual);
%docstring(Domain::allocateJacobian);
%docstring(Domain::allocateSolution);

%docstring(Domain::assemble);

%docstring(Domain::buildState);
%docstring(Domain::buildSolution);



PF4PY_NAMESPACE_END
