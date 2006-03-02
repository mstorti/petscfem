// -*- c++ -*-
// $Id: Error.i,v 1.1.2.2 2006/03/02 21:37:12 rodrigop Exp $

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
    SWIG_exception(SWIG_UnknownError, "unknown exception");
  }
}
