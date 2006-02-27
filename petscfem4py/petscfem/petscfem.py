# This file was created automatically by SWIG 1.3.27.
# Don't modify this file, modify the SWIG interface instead.

import _petscfem

# This file is compatible with both classic and new-style classes.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name) or (name == "thisown"):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


Int    = _petscfem.Int
Float  = _petscfem.Float

class Nodedata(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Nodedata, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Nodedata, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ PyPF::Nodedata instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Nodedata, 'this', _petscfem.new_Nodedata(*args))
        _swig_setattr(self, Nodedata, 'thisown', 1)
    def __del__(self, destroy=_petscfem.delete_Nodedata):
        try:
            if self.thisown: destroy(self)
        except: pass

    def getOption(*args): return _petscfem.Nodedata_getOption(*args)
    def setOption(*args): return _petscfem.Nodedata_setOption(*args)
    def getSize(*args): return _petscfem.Nodedata_getSize(*args)
    def getDim(*args): return _petscfem.Nodedata_getDim(*args)
    def getArray(*args): return _petscfem.Nodedata_getArray(*args)
    def setArray(*args): return _petscfem.Nodedata_setArray(*args)

class NodedataPtr(Nodedata):
    def __init__(self, this):
        _swig_setattr(self, Nodedata, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Nodedata, 'thisown', 0)
        self.__class__ = Nodedata
_petscfem.Nodedata_swigregister(NodedataPtr)

class Elemset(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Elemset, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Elemset, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ PyPF::Elemset instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Elemset, 'this', _petscfem.new_Elemset(*args))
        _swig_setattr(self, Elemset, 'thisown', 1)
    def __del__(self, destroy=_petscfem.delete_Elemset):
        try:
            if self.thisown: destroy(self)
        except: pass

    def getOption(*args): return _petscfem.Elemset_getOption(*args)
    def setOption(*args): return _petscfem.Elemset_setOption(*args)
    def setUp(*args): return _petscfem.Elemset_setUp(*args)
    def getType(*args): return _petscfem.Elemset_getType(*args)
    def getName(*args): return _petscfem.Elemset_getName(*args)
    def getConnectivity(*args): return _petscfem.Elemset_getConnectivity(*args)
    def setConnectivity(*args): return _petscfem.Elemset_setConnectivity(*args)
    def getSize(*args): return _petscfem.Elemset_getSize(*args)
    def setNDof(*args): return _petscfem.Elemset_setNDof(*args)
    def getNDof(*args): return _petscfem.Elemset_getNDof(*args)
    def view(*args): return _petscfem.Elemset_view(*args)

class ElemsetPtr(Elemset):
    def __init__(self, this):
        _swig_setattr(self, Elemset, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Elemset, 'thisown', 0)
        self.__class__ = Elemset
_petscfem.Elemset_swigregister(ElemsetPtr)

class Mesh(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Mesh, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Mesh, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ PyPF::Mesh instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Mesh, 'this', _petscfem.new_Mesh(*args))
        _swig_setattr(self, Mesh, 'thisown', 1)
    def __del__(self, destroy=_petscfem.delete_Mesh):
        try:
            if self.thisown: destroy(self)
        except: pass

    def getOption(*args): return _petscfem.Mesh_getOption(*args)
    def setOption(*args): return _petscfem.Mesh_setOption(*args)
    def getNodedata(*args): return _petscfem.Mesh_getNodedata(*args)
    def setNodedata(*args): return _petscfem.Mesh_setNodedata(*args)
    def addElemset(*args): return _petscfem.Mesh_addElemset(*args)
    def getElemset(*args): return _petscfem.Mesh_getElemset(*args)
    def getSize(*args): return _petscfem.Mesh_getSize(*args)
    def hasElemset(*args): return _petscfem.Mesh_hasElemset(*args)
    def findElemset(*args): return _petscfem.Mesh_findElemset(*args)

class MeshPtr(Mesh):
    def __init__(self, this):
        _swig_setattr(self, Mesh, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Mesh, 'thisown', 0)
        self.__class__ = Mesh
_petscfem.Mesh_swigregister(MeshPtr)

class Constraint(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Constraint, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Constraint, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ PyPF::Constraint instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Constraint, 'this', _petscfem.new_Constraint(*args))
        _swig_setattr(self, Constraint, 'thisown', 1)
    def __del__(self, destroy=_petscfem.delete_Constraint):
        try:
            if self.thisown: destroy(self)
        except: pass


class ConstraintPtr(Constraint):
    def __init__(self, this):
        _swig_setattr(self, Constraint, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Constraint, 'thisown', 0)
        self.__class__ = Constraint
_petscfem.Constraint_swigregister(ConstraintPtr)

class DofMap(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DofMap, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DofMap, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ PyPF::DofMap instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, DofMap, 'this', _petscfem.new_DofMap(*args))
        _swig_setattr(self, DofMap, 'thisown', 1)
    def __del__(self, destroy=_petscfem.delete_DofMap):
        try:
            if self.thisown: destroy(self)
        except: pass

    def addFixations(*args): return _petscfem.DofMap_addFixations(*args)
    def addConstraints(*args): return _petscfem.DofMap_addConstraints(*args)
    def view(*args): return _petscfem.DofMap_view(*args)

class DofMapPtr(DofMap):
    def __init__(self, this):
        _swig_setattr(self, DofMap, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, DofMap, 'thisown', 0)
        self.__class__ = DofMap
_petscfem.DofMap_swigregister(DofMapPtr)

class Fixation(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Fixation, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Fixation, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ PyPF::Fixation instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["node"] = _petscfem.Fixation_node_set
    __swig_getmethods__["node"] = _petscfem.Fixation_node_get
    if _newclass:node = property(_petscfem.Fixation_node_get, _petscfem.Fixation_node_set)
    __swig_setmethods__["field"] = _petscfem.Fixation_field_set
    __swig_getmethods__["field"] = _petscfem.Fixation_field_get
    if _newclass:field = property(_petscfem.Fixation_field_get, _petscfem.Fixation_field_set)
    def __init__(self, *args):
        _swig_setattr(self, Fixation, 'this', _petscfem.new_Fixation(*args))
        _swig_setattr(self, Fixation, 'thisown', 1)
    def __del__(self, destroy=_petscfem.delete_Fixation):
        try:
            if self.thisown: destroy(self)
        except: pass


class FixationPtr(Fixation):
    def __init__(self, this):
        _swig_setattr(self, Fixation, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Fixation, 'thisown', 0)
        self.__class__ = Fixation
_petscfem.Fixation_swigregister(FixationPtr)



