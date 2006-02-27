// -*- c++ -*-

%module petscfem


%include typemaps.i
%include array.i
%rename(view) print;
%include stl.i

%import  namespace.h

%include init.i

%include Error.i
%include SmartPtr.i
%include Nodedata.i
%include Elemset.i
%include Mesh.i
%include Constraint.i
%include DofMap.i
%include Fixation.i
