// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: pfgtpmacr.h,v 1.1 2004/10/24 16:25:21 mstorti Exp $
#ifndef PETSCFEM_PFGTPMACR_H
#define PETSCFEM_PFGTPMACR_H

#define TGETOPTDEF_ND_PFMAT(thash,type,name,default)		\
        name = default;						\
        ierr = ::get_##type(thash,#name,&name,1);		\
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

#define TGETOPTDEF_S_ND_PFMAT(thash,type,name,default)		\
        name=type(#default);					\
        ierr = ::get_##type(thash,#name,name,1);		\
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

#define TGETOPTDEF_S_PFMAT(thash,type,name,default)		\
        type name=type(#default);				\
        ierr = get_##type(thash,#name,name,1);			\
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

#define TGETOPTDEF_ND_PF(thash,type,name,default)		\
        name = default;						\
        get_option(#name,&name); 
  
#define TGETOPTDEF_S_ND_PF(thash,type,name,default)	\
        name = string(#default);			\
        get_option(#name,name); 
  
#define TGETOPTDEF_S_PF(thash,type,name,default)	\
        string name;					\
        TGETOPTDEF_S_ND_PF(thash,type,name,default)
  
#endif
