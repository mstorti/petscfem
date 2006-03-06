// -*- c++ -*-
// $Id: macros.h,v 1.1.2.1 2006/03/06 16:56:04 rodrigop Exp $

#ifndef PYPF_MACROS_H
#define PYPF_MACROS_H


#define PYPF_CLASS(NAME) \
class NAME : public SmartPtr< ::NAME >, public Object


#define PYPF_CTOR_FROM_PTR(NAME) \
public: \
NAME(::NAME* p) : Ptr(p) \
{ if (p==NULL) throw Error("null pointer to "#NAME); } \
private:


#define PYPF_OBJ_GETOPTTBL_DECL \
protected: \
OptionTable* get_opt_table(bool create=false); \
private:


#endif // PYPF_MACROS_H
