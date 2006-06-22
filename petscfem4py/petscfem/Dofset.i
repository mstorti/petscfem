// -*- c++ -*-

%include Object.i
%include Amplitude.i

%apply (int*, int*) { (int* nnod, int* ndof) };

%typemap(doc,name="nodes", type="int const []")    (int nn, const int    nodes[])  "$2_name: $2_type value";
%typemap(doc,name="nodes1",type="int const []")    (int n1, const int    nodes1[]) "$2_name: $2_type value";
%typemap(doc,name="nodes2",type="int const []")    (int n2, const int    nodes2[]) "$2_name: $2_type value";
%typemap(doc,name="fields",type="int const []")    (int nf, const int    fields[]) "$2_name: $2_type value";
%typemap(doc,name="values",type="double const []") (int nv, const double values[]) "$2_name: $2_type value";
%typemap(doc,name="coeffs",type="double const []") (int nc, const double coeffs[]) "$2_name: $2_type value";
ARRAY_FLAT(int nn, const int    nodes[],  ARRAY_INPUT, PyPF_INT)
ARRAY_FLAT(int n1, const int    nodes1[], ARRAY_INPUT, PyPF_INT)
ARRAY_FLAT(int n2, const int    nodes2[], ARRAY_INPUT, PyPF_INT)
ARRAY_FLAT(int nf, const int    fields[], ARRAY_INPUT, PyPF_INT)
ARRAY_FLAT(int nv, const double values[], ARRAY_INPUT, PyPF_FLOAT)
ARRAY_FLAT(int nc, const double coeffs[], ARRAY_INPUT, PyPF_FLOAT)



PYPF_NAMESPACE_BEGIN
%extend Dofset {
  void setFixation(int nn, const int    nodes[],
		   int nf, const int    fields[],
		   int nv, const double values[]) {
    using namespace PYPF_NAMESPACE;
    std::vector<int>    _f;
    std::vector<double> _v;
    if (nn == 0 || nf == 0) return;
    if (nn > 1 && nf == 1) { _f.resize(nn, fields[0]); nf = nn; fields = &_f[0]; }
    if (nn > 1 && nv == 1) { _v.resize(nn, values[0]); nv = nn; values = &_v[0]; }
    PYPF_ASSERT(nn == nf, "incompatible array sizes, size(nodes) != size(fields)");
    PYPF_ASSERT(nn == nv, "incompatible array sizes, size(nodes) != size(values)");
    self->addFixations(nn, nodes, fields, values);
  }
  void setFixation(int nn, const int nodes[],
		   int nf, const int fields[],
		   Amplitude& amplitude) {
    using namespace PYPF_NAMESPACE;
    std::vector<int>    _f;
    std::vector<double> _v;
    if (nn == 0 || nf == 0) return;
    if (nn > 1 && nf == 1) { _f.resize(nn, fields[0]); nf = nn; fields = &_f[0]; }
    PYPF_ASSERT(nn == nf, "incompatible array sizes, size(nodes) != size(fields)");
    std::vector<double> values(nn, 1.0);
    self->addFixations(nn, nodes, fields, &values[0], &amplitude);
  }
  void setFixation(int nn, const int    nodes[],
		   int nf, const int    fields[],
		   int nv, const double values[],
		   Amplitude& amplitude) {
    using namespace PYPF_NAMESPACE;
    std::vector<int>    _f;
    std::vector<double> _v;
    if (nn == 0 || nf == 0) return;
    if (nn > 1 && nf == 1) { _f.resize(nn, fields[0]); nf = nn; fields = &_f[0]; }
    if (nn > 1 && nv == 1) { _v.resize(nn, values[0]); nv = nn; values = &_v[0]; }
    PYPF_ASSERT(nn == nf, "incompatible array sizes, size(nodes) != size(fields)");
    PYPF_ASSERT(nn == nv, "incompatible array sizes, size(nodes) != size(values)");
    self->addFixations(nn, nodes, fields, values, &amplitude);
  }

  void setPeriodic(int n1, const int nodes1[],
		   int n2, const int nodes2[]) {
    using namespace PYPF_NAMESPACE;
    double C[2] = {1.0, -1.0};
    int    N[2];
    int    F[2];
    if ((n1 == 0 && n2 == 0)) return;
    PYPF_ASSERT(n1 == n2, "incompatible array sizes, size(nodes1) != size(nodes2)");
    int ndof; self->getSizes(NULL, &ndof);
    for(int n=0; n<n1; n++) {
      N[0] = nodes1[n];
      N[1] = nodes2[n];
      for(int f=0; f<ndof; f++) {
	F[0] = F[1] = f;
	self->addConstraints(2, C, N, F);
      }
    }
  }
  void setPeriodic(int n1, const int nodes1[],
		   int n2, const int nodes2[],
		   int nf, const int fields[]) {
    using namespace PYPF_NAMESPACE;
    double C[2] = {1.0, -1.0};
    int    N[2];
    int    F[2];
    if ((n1 == 0 && n2 == 0) || nf == 0) return;
    PYPF_ASSERT(n1 == n2, "incompatible array sizes, size(nodes1) != size(nodes2)");
    for(int n=0; n<n1; n++) {
      N[0] = nodes1[n]; 
      N[1] = nodes2[n];
      for(int f=0; f<nf; f++) {
	F[0] = F[1] = fields[f]; 
	self->addConstraints(2, C, N, F);
      }
    }
  }
  void setConstraint(int nc, const double coeffs[],
		     int nn, const int    nodes[],
		     int nf, const int    fields[]) {
    using namespace PYPF_NAMESPACE;
    if (nn == 0 && nf == 0 && nc == 0) return;
    PYPF_ASSERT(nc > 1,   "array size too small, size(coeffs) == 1");
    PYPF_ASSERT(nn > 1,   "array size too small, size(nodes) == 1");
    PYPF_ASSERT(nf > 1,   "array size too small, size(fields) == 1");
    PYPF_ASSERT(nn == nf, "incompatible array sizes, size(nodes) != size(fields)");
    PYPF_ASSERT(nn == nc, "incompatible array sizes, size(nodes) != size(coeffs)");
    self->addConstraints(nc, coeffs, nodes, fields);
  }
}
%ignore Dofset::addFixations;
%ignore Dofset::addConstraints;
PYPF_NAMESPACE_END


%include "Dofset.h"

%clear (int* nnod, int* ndof);

%clear (int nn, const int    nodes[]);
%clear (int n1, const int    nodes1[]);
%clear (int n2, const int    nodes2[]);
%clear (int nf, const int    fields[]);
%clear (int nv, const double values[]);
%clear (int nc, const double coeffs[]);


