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
type name##_12 = default,type name##_13 = default,type name##_14 = default,  \
type name##_15 = default,type name##_16 = default,type name##_17 = default,  \
type name##_18 = default,type name##_19 = default,type name##_20 = default,  \
type name##_21 = default,type name##_22 = default,type name##_23 = default,  \
type name##_24 = default,type name##_25 = default,type name##_26 = default,  \
type name##_27 = default,type name##_28 = default,type name##_29 = default

/** Same, but without the default value
 */
#define ARG_LIST_ND(type,name) \
type name##_0,type name##_1,type name##_2,  \
type name##_3,type name##_4,type name##_5,  \
type name##_6,type name##_7,type name##_8,  \
type name##_9,type name##_10,type name##_11,  \
type name##_12,type name##_13,type name##_14,  \
type name##_15,type name##_16,type name##_17,  \
type name##_18,type name##_19,type name##_20,  \
type name##_21,type name##_22,type name##_23,  \
type name##_24,type name##_25,type name##_26,  \
type name##_27,type name##_28,type name##_29

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
if (name##_15 == default) goto exit_label; \
    indx.push_back(name##_15);          \
                                          \
if (name##_16 == default) goto exit_label; \
    indx.push_back(name##_16);          \
                                          \
if (name##_17 == default) goto exit_label; \
    indx.push_back(name##_17);          \
                                          \
if (name##_18 == default) goto exit_label; \
    indx.push_back(name##_18);          \
                                          \
if (name##_19 == default) goto exit_label; \
    indx.push_back(name##_19);          \
                                          \
if (name##_20 == default) goto exit_label; \
    indx.push_back(name##_20);          \
                                          \
if (name##_21 == default) goto exit_label; \
    indx.push_back(name##_21);          \
                                          \
if (name##_22 == default) goto exit_label; \
    indx.push_back(name##_22);          \
                                          \
if (name##_23 == default) goto exit_label; \
    indx.push_back(name##_23);          \
                                          \
if (name##_24 == default) goto exit_label; \
    indx.push_back(name##_24);          \
                                          \
if (name##_25 == default) goto exit_label; \
    indx.push_back(name##_25);          \
                                          \
if (name##_26 == default) goto exit_label; \
    indx.push_back(name##_26);          \
                                          \
if (name##_27 == default) goto exit_label; \
    indx.push_back(name##_27);          \
                                          \
if (name##_28 == default) goto exit_label; \
    indx.push_back(name##_28);          \
                                          \
if (name##_29 == default) goto exit_label; \
    indx.push_back(name##_29);          \
                                          \
exit_label:;

// This header file automatically generated from readlist.eperl
//                       DON'T MODIFY!!
// In any case, modify the eperl script and re-run eperl
// with make `readlist.h'
