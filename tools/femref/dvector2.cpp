// $Id: dvector2.cpp,v 1.3 2005/01/17 02:47:47 mstorti Exp $
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <libguile.h>

#include <src/dvector.h>

void scmlist2vec(SCM s_list,vector<int> &v) {
  if (!scm_list_p(s_list))
    scm_misc_error ("scmlist2vec",
		    "Not received list.",SCM_EOL);
  v.clear();
  SCM q = s_list;
  while(!SCM_NULLP(q)) {
    SCM item = SCM_CAR(q);
    if (!scm_integer_p(item))
      scm_misc_error ("scmlist2vec",
		      "Not integer value.",SCM_EOL);
    int w = SCM_INUM(item);
    v.push_back(w);
    q = SCM_CDR(q);
  }
}
