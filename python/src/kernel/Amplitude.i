// -*- c++ -*-
// $Id$

PF4PY_NAMESPACE_BEGIN
%director Amplitude;
%director Amplitude::operator();
PF4PY_NAMESPACE_END

PF4PY_NAMESPACE_BEGIN
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
PF4PY_NAMESPACE_END

%include "Amplitude.h"
