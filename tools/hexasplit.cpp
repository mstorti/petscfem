//__INSERT_LICENSE__
//$Id: hexasplit.cpp,v 1.3 2002/07/13 23:19:21 mstorti Exp $
#include <stdio.h>
#include <unistd.h>
#include <vector>
#include <set>
#include <deque>
#include <string>

#define NEL (8) // number of nodes per element
#define NFACES (6) // number of faces per element
#define NELFACE (4) // number of nodes per faces 
#define NSUBEL (5) // number of tetras in an hexa
#define NODSUBEL (4) // number of nodes in each subel

struct row {
  int row_[NEL];
};

enum ElemState {NOT_SET = 0x00000, 
		SPLIT_P = 0x00001, 
		SPLIT_M = 0x00002, 
		SPLIT   = SPLIT_M | SPLIT_P,
		MULTI_SPLIT = 0x00004,
		ENQUEUED = 0x00008};

void e2faces(int k,int *icone_row, set<int> *faces, int nfaces[][NELFACE]) {

  for (int j=0; j<NFACES; j++) faces[j].erase(faces[j].begin(),faces[j].end());

#define INSERT_NODE(face,node) \
 {faces[face].insert(icone_row[node]); \
  nfaces[face][kk++] = icone_row[node];}

  int kk=0;
  INSERT_NODE(0,0);
  INSERT_NODE(0,3);
  INSERT_NODE(0,2);
  INSERT_NODE(0,1);

  kk=0;
  INSERT_NODE(1,0);
  INSERT_NODE(1,1);
  INSERT_NODE(1,5);
  INSERT_NODE(1,4);

  kk=0;
  INSERT_NODE(2,2);
  INSERT_NODE(2,6);
  INSERT_NODE(2,5);
  INSERT_NODE(2,1);

  kk=0;
  INSERT_NODE(3,2);
  INSERT_NODE(3,3);
  INSERT_NODE(3,7);
  INSERT_NODE(3,6);

  kk=0;
  INSERT_NODE(4,0);
  INSERT_NODE(4,4);
  INSERT_NODE(4,7);
  INSERT_NODE(4,3);

  kk=0;
  INSERT_NODE(5,5);
  INSERT_NODE(5,6);
  INSERT_NODE(5,7);
  INSERT_NODE(5,4);

}

void print_face(set<int> face) {
  set<int>::iterator j;
  for (j=face.begin(); j!=face.end(); j++) {
    printf("%d ",*j);
  }
  printf("\n");
}

int main (int argc, char **argv) {
  char c;
  string icone_file = "icone";
  string icone_tetra = "icone_tetra";
  while ((c = getopt (argc, argv, "i:o:")) != -1) {
    switch (c) {
    case 'i':
      icone_file = string(optarg);
      break;
    case 'o':
      icone_tetra = string(optarg);
      break;
    default:
      abort ();
    }
  }
  vector<row> icone;
  // int row[NEL]={0,0,0,0,0,0,0,0};
  // vector<double[3]> xnod;
  
  FILE *ico_file;
  ico_file = fopen(icone_file.c_str(),"r");

  int val;
  int iele=0;
  int nnod=0;
  while (1) {
    int nread = fscanf(ico_file,"%d",&val);
    if (val>nnod) nnod=val;
    if (nread==EOF) break;
    icone.resize(icone.size()+1);
    icone[iele].row_[0]=val;
    for (int k=1; k<=NEL-1; k++) {
      nread = fscanf(ico_file,"%d",&val);
      icone[iele].row_[k] = val;
      if (val>nnod) nnod=val;
    }
    iele++;
  }
  fclose(ico_file);
  int nelem=iele;

  // formar la lista de elementos conectados a un nodo
  vector<set<int> > nod2ele(nnod);
  for (int k=1; k<=nelem; k++) {
    for (int j=0; j<NEL; j++) {
      nod2ele[icone[k-1].row_[j]-1].insert(k);
    }
  }

#if 0
  for (int j=0; j<nnod; j++) {
    printf("nodo: %d -> elems: ",j+1);
    set<int>::iterator jj;
    for (jj=nod2ele[j].begin(); jj!=nod2ele[j].end(); jj++) {
      printf("%d ",*jj);
    }
    printf("\n");
  }
#endif

  // cola de elementos a procesar
  deque<int> eque;
  set<int> *faces = new set<int>[NFACES];
  set<int> *faces_n = new set<int>[NFACES];
  int nfaces[NFACES][NELFACE];
  int nfaces_n[NFACES][NELFACE];
  int *pflag = new int[nelem]; // flags whether the element has been partitioned
  for (int j=0; j<nelem; j++) {
    pflag[j]=NOT_SET;
  }

  int split_p=0,split_m=0,multi=0;
  eque.push_back(1);
  for (int elem_count=0; 1; elem_count++) {
    if (eque.size()==0) {
      // Search if there is an element not partitioned
      for (int j=0; j<nelem; j++) {
	if (pflag[j] == NOT_SET) {
	  printf("Warning: disconnected mesh?? Coninuing with element %d\n",j+1);
	  eque.push_back(j+1);
	  break;
	}
      }
    }

    // There are no more elements to partition
    if (eque.size()==0) goto EXIT;

    // Get next element to partition
    int ele = eque.front();
    eque.pop_front();
    set<int> eneighbors;
    eneighbors.insert(ele);

    int indsplit[NFACES]; // split inducido por loe elementos adyacentes
    for (int j=0; j<NFACES; j++) indsplit[j]=NOT_SET;
    
    e2faces(ele,icone[ele-1].row_,faces,nfaces);
    for (int kk=0; kk<NEL; kk++) {
      int node = icone[ele-1].row_[kk];
      set<int>::iterator jj;
      for (jj=nod2ele[node-1].begin(); jj!=nod2ele[node-1].end(); jj++) {
	int ele_n = *jj; // neighbor element

	if (eneighbors.find(ele_n) != eneighbors.end()) continue; // el elemento ya fue chequeado
	eneighbors.insert(ele_n);
	e2faces(ele_n,icone[ele_n-1].row_,faces_n,nfaces_n);
	
	for (int jface=0; jface<NFACES; jface++) {
	  // printf("Element %d, face %d, nodes ",ele,jface);
	  // print_face(faces[jface]);
	  for (int kface=0; kface<NFACES; kface++) {
	    // printf("       ->   Element %d, face %d, nodes ",ele_n,kface);
	    // print_face(faces_n[kface]);
	    if (faces[jface] == faces_n[kface]) {
	      // printf("matching faces found!\n");
	      if (pflag[ele_n-1] == NOT_SET) {
		// pone en cola
		eque.push_back(ele_n);
		pflag[ele_n-1] = ENQUEUED;
	      } else if ( pflag[ele_n-1] & SPLIT) {
		// Induce splitting en este
		int jjj;
		for (jjj=0; jjj<NELFACE; jjj++) {
		  if (nfaces_n[kface][jjj]==nfaces[jface][0]) break;
		}
		indsplit[jface] = (jjj % 2==0? pflag[ele_n-1] : pflag[ele_n-1] ^ SPLIT);
	      }
	    }
	  }
	}
      }
    }

    int split=NOT_SET;
    // printf("Elemento %d, \n",ele);
    for (int jface=0; jface<NFACES; jface++) {
#if 0
      printf("face %d,  nodes ",jface+1);
      for (int kk=0; kk<NELFACE; kk++) {
	printf("%d ",nfaces[jface][kk]);
      }
      printf(",  split %d\n",indsplit[jface]);
#endif
      if (split!=0 && indsplit[jface]!=0 && split!=indsplit[jface]) {
	printf("Incompatible split! element %d",ele);
	split=MULTI_SPLIT;
      } else if (indsplit[jface] & SPLIT) {
	split = indsplit[jface];
      }
    }

    // If split not induced, set to SPLIT_P
    if (split==NOT_SET) split=SPLIT_P;
    pflag[ele-1]=split;

    switch(split) {
    case SPLIT_P:
      split_p++;
      break;
    case SPLIT_M:
      split_m++;
      break;
    case MULTI_SPLIT:
      multi++;
      break;
    default:
      printf("warning: unknown element splitting %d\n",split);
    }

#if 0
    printf("Elements in queue: ");
    deque<int>::iterator jj;
    for (jj=eque.begin(); jj!=eque.end(); jj++) {
      printf("%d ",*jj);
    }
    printf("\n");
#endif
  }
 EXIT:;
  printf("Statistics: %d split_p, %d split_m, %d multi\nWriting tetra mesh...\n",
	 split_p,split_m,multi);
  fflush(stdout);
  assert(multi==0);

  int tet_split_p[5][4] = { 
  1, 6, 8, 5,
  1, 3, 6, 2,
  3, 8, 6, 7,
  1, 8, 3, 4,
  1, 3, 8, 6
  };

  int tet_split_m[5][4] = {
  2, 7, 5, 6,
  2, 4, 7, 3,
  4, 5, 7, 8,
  2, 5, 4, 1,
  2, 4, 5, 7
  };

  int (*tet_split)[5][4] ;

  FILE *tetras = fopen(icone_tetra.c_str(),"w");
  for (int j=0; j<nelem; j++) {
    tet_split = (pflag[j]==SPLIT_P ? &tet_split_p : &tet_split_m);
    for (int jj=0; jj<NSUBEL; jj++) {
      for (int jjj=0; jjj<NODSUBEL; jjj++) {
	fprintf(tetras,"%d ",icone[j].row_[(*tet_split)[jj][jjj]-1]);
      }
      fprintf(tetras,"\n");
    }
  }
  fclose(tetras);
  printf("Done\n",
}      

