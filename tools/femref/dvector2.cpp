// $Id: dvector2.cpp,v 1.2 2005/01/16 23:59:32 mstorti Exp $
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <libguile.h>

#include <src/dvector.h>

void scmlist2vec(SCM s_list,vector<int> &v) {
  SCM_ASSERT(scm_list_p(s_list),
	     s_list, SCM_ARG1, __FUN__);
  v.clear();
  SCM q = s_list;
  while(!SCM_NULLP(q)) {
    SCM item = SCM_CAR(q);
    SCM_ASSERT(scm_integer_p(item),item,SCM_ARG1,"scmlist2vec");
    int w = SCM_INUM(item);
    printf("%d ",w);
    v.push_back(w);
    q = SCM_CDR(q);
  }
  printf("\n");
}
