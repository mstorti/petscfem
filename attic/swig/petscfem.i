%module petscfem
%{
#include <sles.h>
%}

%include perlmain.i

extern void petscfem_initialize(void);
extern int make_rhs_vector(Vec *);
extern void print_vector(Vec *);
