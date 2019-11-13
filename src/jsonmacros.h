// -*- mode: c++ -*-
#define GETD(X) X = opts.get(#X,NAN).asDouble(); \
  PETSCFEM_ASSERT0(!ISNAN(X),#X " is required")
#define GETI(X) X = opts.get(#X,-1).asInt(); \
  PETSCFEM_ASSERT0(X!=-1,#X " is required") 
// DEFINE version: Declares obtain and check the value
#define GETDD(X) cs_real_t GETD(X)
#define GETID(X) int GETI(X)
// PLAIN version: Just get the value, with a default
#define GETDP(X,DEF) X = opts.get(#X,DEF).asDouble()
#define GETIP(X,DEF) X = opts.get(#X,DEF).asInt()
