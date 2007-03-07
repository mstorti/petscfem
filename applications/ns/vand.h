// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: vand.h,v 1.1.2.1 2007/03/07 00:56:14 mstorti Exp $
#ifndef PETSCFEM_VAND_H
#define PETSCFEM_VAND_H

extern FILE* dump_file;
extern dvector<int> vd_elems_loc,vd_ebuff,vd_elems;
extern dvector<double> vd_data_loc,vd_dbuff,vd_data;
#define VD_DATA_SIZE 4

#endif
