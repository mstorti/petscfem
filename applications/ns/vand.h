// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: vand.h,v 1.1.4.3 2007/03/15 02:49:11 mstorti Exp $
#ifndef PETSCFEM_VAND_H
#define PETSCFEM_VAND_H

#define VD_DATA_SIZE 4
#include <map>

struct VDDumpData {
  dvector<int> vd_elems_loc;
  dvector<double> vd_data_loc;
};

typedef map <string,VDDumpData *> vd_map_t;
extern vd_map_t vd_map;
extern int vd_dump_flag;

#endif
