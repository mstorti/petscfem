#include <glib.h>

#include "fem.h"
#include "nodedata.h"

#define ML 20
void int2s(int n,string &s) {
  char nmbr[ML]; // This is ugly!!
  snprintf(nmbr,ML,"%d",j+1);
  s=string(nmbr);
}

void build_map(SList &sl, SMap &smap) {
  smap.clear();
  int jj=0;
  for (SList::iterator j=sl.begin(); j!=sl.end(); j++) {
    smap[*j]=jj++;
  }
}

NodeData::NodeData(FileStack &fstack) {
  part=1;

  char *line;
  const char *cline;
  read_hash_table(fstack,thash);
  
  thash->get_entry("fields",cline);
  if (cline) {
    parse_fields_line(field_list,line);
    thash->get_entry("equations",cline);
    if(cline) {
      parse_fields_line(eq_list,line);
      build_map(eq_list,eq_map);
    } else {
      eq_list = field_list;
      eq_map = fld_map;
    }
  } else {
    int ndof=0;
    get_int(thash,"ndof",ndof,1);
    PETSCFEM_ERROR(ndof==0,"Not \"ndof\" or \"fields\" line in texthash \n"
		   "table for \"nodes\"");
    string s;
    for (int j=1; j<=ndof, j++) {
      int2s(j,s);
      s.insert(0,"u[").append("]");
      field_list.push_back(s);
    }
    build_map(field_list,field_map);
    eq_list = field_list;
    eq_map = fld_map;
  }

  ndof = field_list.size();
  PETSCFEM_ERROR(ndof<=0,"\"ndof\" not positive while "
		 "constructing NodeData\n");
  neqs = eq_list.size();

  nnod = 0;
  get_int(thash,"nnod",nnod,1);
  darray *xnod;
  GArray *gdata = g_array_new(0,0,ndof*sizeof(double));
  if (nnod>0) gdata = g_array_set_size(gdata,nnod);
  int nnod_read=nnod;

  char *token;
  int node=0;
  while (1) {
    fstack->get_line(line);
    if (strstr("__END_NODES__",line)) break;

    node++;
    if (node>nnod) {
      nnod=node;
      gdata = g_array_set_size(gdata,nnod);
    }
    double *row = &g_array_index(gdata,double,node);
    for (int kk=0; kk<nu; kk++) {
      token = strtok((kk==0 ? line : NULL),bsp);
      PETSCFEM_ERROR(token==NULL,
		     "Error reading coordinates in line:\n\"%s\"\n"
		     "Not enough values in line!!\n",line);
      int nread = sscanf(token,"%lf",row+kk);
      PETSCFEM_ERROR(nread != 1,
		     "Error reading coordinates in line:\n\"%s\"",line);
    }
  }
  
  PETSCFEM_ERROR(nnod_read>0 && nnod!=nnod_read,
		 "NodeData: error: \"nnod = %d\" entered but "
		 "%d nodes read\n",nnod_read,nnod);

  data = &g_array_index(gdata,double,1);
  PetscPrintf(PETSC_COMM_WORLD,"Read %d nodes\n",nnod);
      
  delete[] row;
  mesh->nodedata->nodedata = new double[nnod*nu];
  for (node=0; node<nnod; node++) {
    row = (double *) da_ref(xnod,node);
    for (int kk=0; kk<nu; kk++) {
      NODEDATA(node,kk) = row[kk];
    }
  }
  da_destroy(xnod);

}

#if 0
void NodeData::set_part_mode(int part_) {
  part=part_;
  if (part) {
    npart = new int[nnod];
  }
}

/// returns the processor where this node (base 0) lives.
int NodeData::proc(int node) {
  if (part) {
    return npart[node];
  } else {
    return 1;
  }
}
#endif
