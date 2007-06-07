// -*- c++ -*-


%template() std::pair<int,int>;
%template() std::vector<int>;
%template() std::vector<float>;

ARRAY_STDVEC_OUTPUT(std::vector<int>& gdofs, PyPF_INT)
ARRAY_STDVEC_OUTPUT(std::vector<int>& ldofs, PyPF_INT)
ARRAY_STDVEC_OUTPUT(std::vector<int>& fdofs, PyPF_INT)

%include "Dofset.h"

%clear std::vector<int>& gdofs;
%clear std::vector<int>& ldofs;
%clear std::vector<int>& fdofs;
