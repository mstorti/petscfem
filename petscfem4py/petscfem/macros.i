// -*- c++ -*-
// $Id: macros.i,v 1.1.2.2 2006/03/06 16:56:04 rodrigop Exp $

#define PYPF_CLASS(CLSNAME) class CLSNAME: public Object

#define PYPF_CTOR_FROM_PTR(CLSNAME)

#define PYPF_OBJ_GETOPTTBL_DECL \
protected: \
OptionTable* get_opt_table(bool create=false); \
private:
