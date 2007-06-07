// -*- c++ -*-
// $Id$

%wrapper %{
#if defined(SWIG_DIRECTORS)
#define PF4PY_CATCH_DIRECTORS\
  catch(Swig::DirectorException &e) { SWIG_fail; }
#else
#define PF4PY_CATCH_DIRECTORS
#endif
#define PF4PY_CATCH_PETSCFEM\
  catch(const PF4PY_NAMESPACE::Error & error)\
  { SWIG_exception(SWIG_RuntimeError, error); } 
#define PF4PY_CATCH_STDEXC\
  catch (std::exception& e)\
  { SWIG_exception(SWIG_SystemError, e.what()); }
#define PF4PY_CATCH_UNKNOWN\
  catch (...)\
  { SWIG_exception(SWIG_UnknownError, "unknown exception"); }
%}

%include exception.i

%exception {
try { $action }
PF4PY_CATCH_DIRECTORS
PF4PY_CATCH_PETSCFEM
PF4PY_CATCH_STDEXC
PF4PY_CATCH_UNKNOWN
}
