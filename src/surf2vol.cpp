//__INSERT_LICENSE__
// $Id: surf2vol.cpp,v 1.5 2003/02/26 01:46:33 mstorti Exp $

#include <src/utils.h>
#include <src/surf2vol.h>
#include <src/surf2vol2.h>
#include <src/linkgraph.h>

extern Mesh *GLOBAL_MESH;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Contains the three non-trivial face orientations. All
    other may be obtained by reflections ad rotations. */
const int Quad2Hexa::faces[][8] = {
  0,1,2,3,4,5,6,7,
  1,5,6,2,0,4,7,3,
  0,4,5,1,3,7,6,2};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Remaps the volume connectivity #vol_conn# so that the face
    is in a standard position. 
    @param surf_map (input) the nodes of the volume element
    that are on the face. 
    @param vol_conn (input/output) the connectivity of the
    volume element. Is remapped so that the face is in a
    standard position. 
*/ 
int Surf2Vol::map_mask(const int *surf_map,int *vol_conn) {
  int nel_surf, nel_vol;
  int nf = nfaces(nel_surf,nel_vol);
  int match=0;
  const int *fc, *vol;
  // Loop over all face orientation posibilities until
  // one of them matches
  for (int f=0; f<nf; f++) {
    face(f,fc,vol);
    match = 1;
    for (int l=0; l<nel_surf; l++) {
      if (surf_map[l] != fc[l]) {
	match=0;
	break;
      }
    }
    if (match) break;
  }
  // Verify that one of the orientations must match
  if (!match) return 0;
  // Remap volume connectivity
  vector<int> vol_conn_c(nel_vol);
  for (int j=0; j<nel_vol; j++) vol_conn_c[j] = vol_conn[vol[j]];
  for (int j=0; j<nel_vol; j++) vol_conn[j] = vol_conn_c[j];
  return 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Quad2Hexa::face(int j,const int *&fc,const int *&vol_ret) { 
  // Changes quad numbering orientation. 
  int spin_map[] = {0,3,2,1};
  // Splits face number j into a triplet: spin(0/1) - rota (0-3) - k (0-2)
  int spin,rota,k,m;
  spin = modulo(j,2,&m);
  rota = modulo(m,4,&k);
  // Construct volume connectivity corresponding to k/rota
  for (int l=0; l<4; l++) {
    int ll = modulo(l+rota,4);
    vol[l] = faces[k][ll];
    vol[l+4] = faces[k][ll+4];
  }
  // Change orientation if needed. 
  if (spin) {
    for (int l=0; l<4; l++) {
      vol_r[l] = vol[4+spin_map[l]];
      vol_r[4+l] = vol[spin_map[l]];
    }
    vol_ret = vol_r;
  } else vol_ret = vol;

  // Face are the first 4 nodes.
  // Rotate depending on `use_exterior_normal'
  for (int k=0; k<4; k++) 
    this_face[k] = (use_exterior_normal() ? 
		    vol_ret[spin_map[k]] : vol_ret[k]);
  fc = this_face;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Line2Quad::face(int j,const int *&fc,const int *&vol) {
  static int fc_c[2], vol_c[4], fc_cr[2], vol_cr[4];
  static const int vol_cc[4] = {0, 1, 3, 2};
  int fc_rot[] = {1,0};
  for (int k=0; k<2; k++) fc_c[k] = (j+k) % 4;
  for (int k=0; k<4; k++) vol_c[k] = (vol_cc[k]+j) % 4;
  if (use_exterior_normal()) {
    fc = fc_c;
    vol = vol_c;
  } else {
    fc_cr[1] = fc_c[0];
    fc_cr[0] = fc_c[1];

    vol_cr[1] = vol_c[0];
    vol_cr[0] = vol_c[1];
    vol_cr[2] = vol_c[3];
    vol_cr[3] = vol_c[2];
    fc = fc_cr;
    vol = vol_cr;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Tri2Prism::rotate(int *vol,int nrot) {
  for (int j=0; j<3; j++) {
    int jj = (j+nrot) % 3;
    vol_c[j] = vol[jj];
    vol_c[3+j] = vol[3+jj];
  }
  for (int j=0; j<6; j++) vol[j] = vol_c[j];
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Tri2Prism::reflect(int *vol) {
  for (int j=0; j<3; j++) {
    int x = vol[3+j];
    vol[3+j] = vol[j];
    vol[j] = x;
  }
  invert(vol);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Tri2Prism::invert(int *vol) {
  int x = vol[0];
  vol[0] = vol[1];
  vol[1] = x;
  x = vol[3];
  vol[3] = vol[4];
  vol[4] = x;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Tri2Prism::face(int j,const int *&fc,const int *&vol_a) {
  int rota, face;
  rota = modulo(j,3,&face);
  for (int j=0; j<6; j++) vol[j] = j;
  rotate(vol,rota);
  if (face) reflect(vol);
  if (use_exterior_normal()) invert(vol);
  fc = vol_a = vol;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Surf2Vol::Surf2Vol(const char *geom,int ndim,int nel,
		   int npg,int mat_version=GP_NEWMAT,
		   int use_exterior_normal_a=0) 
  : GPdata(geom,ndim,nel,npg,mat_version),
    use_exterior_normal_m(use_exterior_normal_a) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void identify_volume_elements_fun(int nnod, int nel_surf, int layers,
			      int nelem, int *icone, int nel, int nel_vol,
			      Elemset *vol_elem,Surf2Vol *sv_gp_data) {
  assert(2*nel_surf==nel_vol);
  for (int layer=0; layer<layers; layer++) {
    // Mark nodes on the surface
    // surface:= is surface[k]==0 then k is not on the surface
    // if != 0 then surface[k] is the number of surface node +1
    vector<int> surface(nnod,0);
    // maps surface numbering (0 to surf_nodes-1) to global (0 to nnod-1)
    vector<int> srf2glb;
    for (int e=0; e<nelem; e++) {
      int *icorow = icone + nel*e + layer* nel_surf;
      for (int j=0; j<nel_surf; j++) surface[icorow[j]-1]=1;
    }
    // Count surface nodes
    int surf_nodes = 0;
    for (int k=0; k<nnod; k++) {
      if (surface[k]) {
	surface[k] = ++surf_nodes;
	srf2glb.push_back(surf_nodes);
      }
    }

    // Construct graph for volume elemset
    LinkGraph graph;
    graph.set_chunk_size(10000);
    graph.init(surf_nodes);

    // Construct node to element array for the volume elemset
    for (int e=0; e<vol_elem->nelem; e++) {
      int *icorow = vol_elem->icone + vol_elem->nel*e;
      for (int j=0; j<nel_vol; j++) {
	int node = icorow[j]-1;
	int snode = surface[node]-1;
	if (snode>=0) graph.add(snode,e);
      }
    }

    // For each surface element look for the corresponding
    // volume element that shares a face
    vector<int> mask(nel_vol), icorow_c(nel_vol);
    for (int e=0; e<nelem; e++) {
      int *icorow = icone + nel*e + layer*nel_surf;
      LinkGraphRow row;
      assert(nel_surf>0);
      // Take list for first node
      int node = icorow[0]-1;
      int sf_node = surface[node]-1;
      graph.set_ngbrs(sf_node,row);
      LinkGraphRow::iterator q;
      int found=0;
      int *vicorow;
      int match=0;
      for (q=row.begin(); q!=row.end(); q++) {
	int ve = *q; // the volume element
	vicorow = vol_elem->icone + vol_elem->nel*ve;
	for (int j=0; j<nel_vol; j++) mask[j]=-1;
	found=0;
	for (int j=0; j<nel_surf; j++) {
	  int sf_node = icorow[j];
	  for (int k=0; k<nel_vol; k++) {
	    if (vicorow[k]==sf_node) {
	      mask[j] = k;
	      break;
	    }
	  }
	  if (mask[j]==-1) break;
	  found++;
	}
	if (found==nel_surf) {
	  for (int j=0; j<nel_vol; j++) icorow_c[j] = vicorow[j];
	  // Volume element was found, find map and
	  if (sv_gp_data->map_mask(mask.begin(),icorow_c.begin())) {
	    match=1;
	    for (int j=0; j<nel_vol; j++) icorow[j] = icorow_c[j];
	    break;
	  }
	}
      }
      if (!match) {
	PetscPrintf(PETSC_COMM_WORLD,
		    "embedded_gatherer: Can't find matching volume element"
		    " to surface element %d\n",e);
	PetscFinalize();
	exit(0);
      }
    }

#if 0
    if (MY_RANK==0) {
      printf("Surface element connectivities: \n");
      for (int e=0; e<nelem; e++) {
	int *icorow = icone + nel*e;
	printf("surf.el. %d: ",e+1);
	for (int j=0; j<nel_vol; j++) printf("%d ",icorow[j]);
	printf("\n");
      }
    }
#endif
 
    graph.clear();
    surface.clear();
    srf2glb.clear();
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Surf2Vol::factory(TextHashTable *thash, string &volume_elemset,
		       int nel_surf,
		       Surf2Vol *&sv_gp_data, Elemset *& vol_elem,
		       int &identify_volume_elements,int &layers,
		       int &use_exterior_normal,int &ndimel) {
  int ierr;
  // Find volume elemset
  vol_elem = GLOBAL_MESH->find(volume_elemset);
  // Verifiy that the elemset was found
  PETSCFEM_ASSERT(vol_elem,"Can't find volume element name: %s\n",
		  volume_elemset.c_str())

  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  //o Number of Gauss points.
  TGETOPTNDEF(thash,int,npg,none);
  // ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  TGETOPTNDEF(thash,int,ndim,none); //nd
  //o Use exterior or interior normal
  TGETOPTDEF_ND(thash,int,use_exterior_normal,1);
  //o Identify automatically the internal volume elements with a face
  // on the surface
  TGETOPTDEF_ND(thash,int,identify_volume_elements,0);
  //o Number of layers in the normal direction.
  TGETOPTDEF_ND(thash,int,layers,2);
  PETSCFEM_ASSERT0(layers>=1,
		   "embedded_gatherer: Number of layers must be integer >=1\n");
  PETSCFEM_ASSERT(layers<=3,"embedded_gatherer: not supported yet layers>2,"
		  " entered layers: %d\n",layers);

  int nel = nel_surf*(layers+1);
  ndimel=ndim-1;
  if (geometry=="quad2hexa") {
    sv_gp_data = new Quad2Hexa(geometry.c_str(),ndim,nel,npg,
			       GP_FASTMAT2,use_exterior_normal);
  } else if (geometry=="tri2prism") {
    sv_gp_data = new Tri2Prism(geometry.c_str(),ndim,nel,npg,
			       GP_FASTMAT2,use_exterior_normal);
  } else if (geometry=="line2quad") {
    sv_gp_data = new Line2Quad(geometry.c_str(),ndim,nel,npg,
			       GP_FASTMAT2,use_exterior_normal);
  } else PETSCFEM_ERROR("embedded_gatherer: unknown geometry \"%s\"\n",geometry.c_str());
}
