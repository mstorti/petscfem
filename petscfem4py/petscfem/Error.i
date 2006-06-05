// -*- c++ -*-
// $Id: Error.i,v 1.1.2.5 2006/06/05 16:02:54 dalcinl Exp $

%include exception.i

%define PYPF_CATCH_SWIGDIRECTORS
/* SWIG Directors exception*/
catch(Swig::DirectorException &e) 
 { SWIG_fail; }
%enddef

%define PYPF_CATCH_PETSCFEM
/* PETScFEM exception */
catch(const PYPF_NAMESPACE::Error & error) 
 { SWIG_exception(SWIG_RuntimeError, error); } 
%enddef

%define PYPF_CATCH_UNKNOWN
/* Other exceptions */
catch (std::exception& e)
 { SWIG_exception(SWIG_SystemError, e.what() ); }
catch (...)
 { SWIG_exception(SWIG_UnknownError, "unknown exception");}
%enddef


%exception {
try { $action }
%#if defined(SWIG_DIRECTORS)
PYPF_CATCH_SWIGDIRECTORS
%#endif
PYPF_CATCH_PETSCFEM
PYPF_CATCH_UNKNOWN
}
