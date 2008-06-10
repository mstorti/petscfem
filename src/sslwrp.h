// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id$
#ifndef PETSCFEM_SSLWRP_H
#define PETSCFEM_SSLWRP_H

// This wrappers must be used in order to avoid warnings
// about string -> char* castings
#define Sopen_wrp(a1,a2) Sopen((char *)(a1),(a2)) 
#define Sprintf_wrp(srvr,s,args...) Sprintf((srvr),(char *)s, ##args) 


#endif
