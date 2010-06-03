# $Id: Makefile,v 1.2 2003/09/02 02:27:37 mstorti Exp $ 

tests: idmap

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Idmap class
idmap: tidmap.bin
#	tidmap.bin >tidmap.sal
	cat /dev/null > tidmap.sal
	./tidmap.bin    5 100 1 >> tidmap.sal
	./tidmap.bin  100  50 0 >> tidmap.sal

tidmap.bin: tidmap.o $(LIBPETSCFEM)
	${CXX_CLINKER} -o $@ tidmap.o $(LDFLAGS_NO_PROG_LIB)

SRCS = 
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

#s
