//__INSERT_LICENSE__
//$Id: hook.cpp,v 1.1 2002/09/23 14:15:11 mstorti Exp $

#include <src/fem.h>
#include <src/readmesh.h>
#include "./nsi_tet.h"

void HookList::init(Mesh &mesh,Dofmap &dofmap) {
  HookList::iterator q;
  for (q=begin(); q!=end(); q++) 
    (*q)->init(mesh,dofmap);
}

void HookList::time_step_pre(double time) {
  HookList::iterator q;
  for (q=begin(); q!=end(); q++) 
    (*q)->time_step_pre(time);
}

void HookList::time_step_post(double time) {
  HookList::iterator q;
  for (q=begin(); q!=end(); q++) 
    (*q)->time_step_post(time);
}
