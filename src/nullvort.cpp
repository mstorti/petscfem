//__INSERT_LICENSE__
// $Id: nullvort.cpp,v 1.2 2003/02/25 20:34:22 mstorti Exp $

#include <src/nullvort.h>
#include <src/dvector.h>
#include <src/surf2vol.h>
#include <src/util2.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void null_vort::read(FileStack *fstack,Mesh *mesh,Dofmap *dofmap) {
  thash.read(fstack);
  int ierr;
  //o The number of nodes per skin panel. 
  TGETOPTNDEF(&thash,int,nel_surf,0);
  //o The elemset that is from the fluid side. 
  TGETOPTDEF_S(&thash,string,volume_elemset,<none>);
  assert(volume_elemset!="<none>");

  int identify_volume_elements, layers,
    use_exterior_normal,  ndimel;
  Surf2Vol *sv_gp_data;
  Elemset *vol_elem;
  Surf2Vol::factory(&thash, volume_elemset, sv_gp_data, 
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

}
