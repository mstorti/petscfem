//__INSERT_LICENSE__
//$Id: hook.cpp,v 1.3 2002/09/23 21:18:16 mstorti Exp $

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include "./nsi_tet.h"

extern int MY_RANK,SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class rosi_hook : public Hook {
private:
  int petscfem2pfm_verbose_m;
  FILE *petscfem2pfm;
public:
  void init(Mesh &mesh,Dofmap &dofmap,const char *name);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void rosi_hook::init(Mesh &mesh,Dofmap &dofmap,const char *name) {
  // File to send forces and moments to "PFM"
  int ierr;
  TextHashTable *t = mesh.global_options;
  TGETOPTDEF_S(t,string,petscfem2pfm_file,);
  TGETOPTDEF(t,int,petscfem2pfm_verbose,0);
  petscfem2pfm=NULL;
  if (petscfem2pfm_file != "" && !MY_RANK) {
    petscfem2pfm = fopen(petscfem2pfm_file.c_str(),"w");
    assert(petscfem2pfm);
    setvbuf(petscfem2pfm,NULL,_IOLBF,0);
  }    
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void rosi_hook::time_step_post(double time,int step,
			       const vector<double> &gather_values) {
  if (petscfem2pfm && !MY_RANK) {
    assert(gather_values.size() >= 6);
    if (petscfem2pfm_verbose_m) {
      printf("rosi_hook: sending to PFM: time %f, forces %f %f %f\n",
	     time,gather_values[0],gather_values[1],gather_values[2]);
      printf("    moments %f %f %f\n",
	     gather_values[3],gather_values[4],gather_values[5]);
    }
    fprintf(petscfem2pfm,"time %e\n",time);
    fprintf(petscfem2pfm,"forces %e %e %e\n",
	    gather_values[0],gather_values[1],gather_values[2]);
    fprintf(petscfem2pfm,"moments %e %e %e\n",
	    gather_values[3],gather_values[4],gather_values[5]);
    fflush(petscfem2pfm);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void HookList::init(Mesh &mesh,Dofmap &dofmap) {
  Hook *hook;
  const char *line;
  mesh.global_options->get_entry("hook_list",line);
  if (!line) return;
  char *lcpy = local_copy(line);
  int n=0; 
  while (1) { 
    char *token = strtok((n++ == 0 ? lcpy : NULL),"[ \t\n]");
    if (!token) break;
    if (!strcmp(token,"rosi_hook")) hook = new rosi_hook;
    else PetscPrintf(PETSC_COMM_WORLD,
		     "Unknown hook \"%s\nLine: \"%s\"\n",
		     token,line);
    PETSCFEM_ASSERT(hook,"Couldn't create hook \"%s\"\n",token);

    token = strtok((n++ == 0 ? lcpy : NULL),"[ \t\n]");
    hook->init(mesh,dofmap,token);
    push_back(hook);
  }
  delete[] lcpy;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void HookList::time_step_pre(double time,int step) {
  HookList::iterator q;
  for (q=begin(); q!=end(); q++) 
    (*q)->time_step_pre(time,step);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void HookList::time_step_post(double time,int step,
			      const vector<double> &gather_values) {
  HookList::iterator q;
  for (q=begin(); q!=end(); q++) 
    (*q)->time_step_post(time,step,gather_values);
}

