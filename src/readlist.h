// -*- mode: c++ -*-
// This header file automatically generated from readlist.eperl
//                       DON'T MODIFY!!
// In any case, modify the eperl script and re-run eperl
// with make `readlist.h'

/** Defines and argument list to be used as a variable argument list.
    @doc Example: ARG\_LIST(int,arg,0) expands to 
    int arg\_0 = 0, int arg\_1 = 0, int arg\_2 = 0 ...
    The current maximum number of args is \$maxargs=.
    This is set in `readlist.eperl'. 
    @author M. Storti
    @param type (input) the type of the variable argument list (int,
    double, etc...)
    @param name (input) the name of the arguments.
    @param default (input) the default value. 
*/ 
#define ARG_LIST(type,name,default) \
type name##_0 = default,type name##_1 = default,type name##_2 = default,  \
type name##_3 = default,type name##_4 = default,type name##_5 = default,  \
type name##_6 = default,type name##_7 = default,type name##_8 = default,  \
type name##_9 = default,type name##_10 = default,type name##_11 = default,  \
type name##_12 = default,type name##_13 = default,type name##_14 = default

/** Same, but without the default value
 */
#define ARG_LIST_ND(type,name) \
type name##_0,type name##_1,type name##_2,  \
type name##_3,type name##_4,type name##_5,  \
type name##_6,type name##_7,type name##_8,  \
type name##_9,type name##_10,type name##_11,  \
type name##_12,type name##_13,type name##_14

/** Reads the argument list in a FastVector of the corresponding type. 
    @doc Example: READ\_ARG\_LIST(name,indx,default,exit\_label)
    This is set in `readlist.eperl'. 
    @author M. Storti
    @param name (input) the name of the arguments.
    @param indx (input) the FastVector where arguments are stored.
    @param default (input) the default value. 
    @param exit\_label (input) the exit label to be generated. Usually: EXIT
*/ 
#define READ_ARG_LIST(name,indx,default,exit_label) \
if (name##_0 == default) goto exit_label; \
    indx.push_back(name##_0);          \
                                          \
if (name##_1 == default) goto exit_label; \
    indx.push_back(name##_1);          \
                                          \
if (name##_2 == default) goto exit_label; \
    indx.push_back(name##_2);          \
                                          \
if (name##_3 == default) goto exit_label; \
    indx.push_back(name##_3);          \
                                          \
if (name##_4 == default) goto exit_label; \
    indx.push_back(name##_4);          \
                                          \
if (name##_5 == default) goto exit_label; \
    indx.push_back(name##_5);          \
                                          \
if (name##_6 == default) goto exit_label; \
    indx.push_back(name##_6);          \
                                          \
if (name##_7 == default) goto exit_label; \
    indx.push_back(name##_7);          \
                                          \
if (name##_8 == default) goto exit_label; \
    indx.push_back(name##_8);          \
                                          \
if (name##_9 == default) goto exit_label; \
    indx.push_back(name##_9);          \
                                          \
if (name##_10 == default) goto exit_label; \
    indx.push_back(name##_10);          \
                                          \
if (name##_11 == default) goto exit_label; \
    indx.push_back(name##_11);          \
                                          \
if (name##_12 == default) goto exit_label; \
    indx.push_back(name##_12);          \
                                          \
if (name##_13 == default) goto exit_label; \
    indx.push_back(name##_13);          \
                                          \
if (name##_14 == default) goto exit_label; \
    indx.push_back(name##_14);          \
                                          \
exit_label:;

// This header file automatically generated from readlist.eperl
//                       DON'T MODIFY!!
// In any case, modify the eperl script and re-run eperl
// with make `readlist.h'
