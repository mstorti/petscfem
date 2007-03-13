// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: vand.h,v 1.1.4.2 2007/03/13 02:50:34 mstorti Exp $
#ifndef PETSCFEM_VAND_H
#define PETSCFEM_VAND_H

#define VD_DATA_SIZE 4
#include <map>

struct VDDumpData {
  int vd_dump_flag;
  dvector<int> vd_elems_loc,vd_ebuff,vd_elems;
  dvector<double> vd_data_loc,vd_dbuff,vd_data;
};

typedef map <string,VDDumpData *> vd_map_t;
extern vd_map_t vd_map;
extern int vd_dump_flag;

#endif
