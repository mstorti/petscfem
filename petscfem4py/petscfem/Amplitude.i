// -*- c++ -*-
// $Id: Amplitude.i,v 1.1.2.5 2006/06/22 22:35:42 dalcinl Exp $

PYPF_NAMESPACE_BEGIN

%feature("ref")   Amplitude "$this->incref();"
%feature("unref") Amplitude "$this->decref();"

%director Amplitude;
%director Amplitude::operator();

PYPF_NAMESPACE_END

PYPF_NAMESPACE_BEGIN
%extend Amplitude {
  %pythoncode {
  def __new__(cls, *args, **kwargs):
      amp = super(Amplitude, cls).__new__(cls)
      if cls != Amplitude:
          kind = getattr(cls, 'KIND', None)
          if kind is None:
              Amplitude.__init__(amp)
          else:
              if isinstance(kind, str):
                  kind = getattr(Amplitude, kind.upper())
              Amplitude.__init__(amp, kind)
      return amp
  }
}
PYPF_NAMESPACE_END

%include "Amplitude.h"
