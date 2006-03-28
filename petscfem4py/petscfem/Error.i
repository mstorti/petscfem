// -*- c++ -*-
// $Id: Error.i,v 1.1.2.3 2006/03/28 22:13:25 rodrigop Exp $

%include exception.i

#if 1

%exception {
  try {
    $action
  } catch(const PyPF::Error & error) {
    SWIG_exception(SWIG_RuntimeError, error);
  } catch(const std::exception& exc) {
    SWIG_exception(SWIG_RuntimeError, exc.what());
  } catch(...) {
    SWIG_exception(SWIG_UnknownError, "unknown exception");
  }
}

#else

%{
#define PF_TRY(ACTION)                             \
try { ACTION }                                     \
catch(const PyPF::Error & error)                   \
{ SWIG_exception(SWIG_RuntimeError, error); }      \
catch(const std::exception& exc)                   \
{ SWIG_exception(SWIG_RuntimeError, exc.what()); } \
catch(...)                                         \
{ SWIG_exception(SWIG_UnknownError, "unknown exception"); }
%}

%exception { PF_TRY($action) }

#endif
