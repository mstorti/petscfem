// -*- c++ -*-
// $Id: Error.i,v 1.1.2.6 2006/06/08 16:05:08 dalcinl Exp $

%wrapper %{
#if defined(SWIG_DIRECTORS)
#define PYPF_CATCH_DIRECTORS\
  catch(Swig::DirectorException &e) { SWIG_fail; }
#else
#define PYPF_CATCH_DIRECTORS
#endif
#define PYPF_CATCH_PETSCFEM\
  catch(const PYPF_NAMESPACE::Error & error)\
  { SWIG_exception(SWIG_RuntimeError, error); } 
#define PYPF_CATCH_UNKNOWN\
  catch (std::exception& e)\
  { SWIG_exception(SWIG_SystemError, e.what() ); }\
  catch (...)\
  { SWIG_exception(SWIG_UnknownError, "unknown exception"); }
%}

%include exception.i

%exception {
try { $action }
PYPF_CATCH_DIRECTORS
PYPF_CATCH_PETSCFEM
PYPF_CATCH_UNKNOWN
}
