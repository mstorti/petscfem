// -*- c++ -*-

%{
#include <exception>
#include "Error.h"
%}

%include exception.i

%exception {
  try {
    $action
  } 
  catch(const PyPF::Error & error) {
    SWIG_exception(SWIG_RuntimeError, error);
  }
  catch(const std::exception& exc) {
    SWIG_exception(SWIG_RuntimeError, exc.what());
  }
  catch(...) {
    SWIG_exception(SWIG_RuntimeError,"Unknown exception");
  }
}
