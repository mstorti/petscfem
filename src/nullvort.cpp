//__INSERT_LICENSE__
// $Id: nullvort.cpp,v 1.4 2003/02/26 01:46:33 mstorti Exp $

#include <src/nullvort.h>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <src/surf2vol.h>
#include <src/util2.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
null_vort::null_vort() { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
null_vort::~null_vort() { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void null_vort::read(FileStack *fstack,Mesh *mesh,Dofmap *dofmap) {
  // Read options from data file
  thash.read(fstack);
  int ierr;

  // get options. Other options are processed in Surf2vol::factory
  //o The number of nodes per skin panel. 
  TGETOPTNDEF(&thash,int,nel_surf,0);
  //o The elemset that is from the fluid side. 
  TGETOPTDEF_S(&thash,string,volume_elemset,<none>);
  assert(volume_elemset!="<none>");

  // Call the `factory' (constructor) for the Surf2Vol object. 
  int identify_volume_elements, layers,
    use_exterior_normal,  ndimel;
  Surf2Vol *sv_gp_data=NULL;
  Elemset *vol_elem;
  Surf2Vol::factory(&thash, volume_elemset, 
		    nel_surf, sv_gp_data, 
		    vol_elem, identify_volume_elements, layers,
		    use_exterior_normal,  ndimel);
  int nel = nel_surf*(layers+1);

  dvector<int> icone;
  icone.reshape(2,0,nel);
  vector<string> tokens;
  vector<int> row;
  int nelem=0;
  char *line;
  while(!fstack->get_line(line)) {
    if (!strcmp(line,"__END_DATA__")) break;
    row.clear();
    read_int_array(row,line);
    assert(row.size()==nel_surf);
    for (int j=0; j<nel_surf; j++) icone.push(row[j]);
    for (int j=0; j<nel_surf*layers; j++) icone.push(1);
    nelem++;
  }
  icone.defrag();
  assert(nelem*nel == icone.size());
  assert(volume_elemset!="<none>");

  identify_volume_elements_fun(dofmap->nnod, nel_surf, layers,
			       nelem, icone.buff(), nel, vol_elem->nel,
			       vol_elem, sv_gp_data);
  delete sv_gp_data;

}
