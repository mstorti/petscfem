PETSCFEM_DIR = ..

SRCS = tryme.cpp getarraypars.tab.c

include $(PETSCFEM_DIR)/Makefile.base

tryme: LINK.o = g++
tryme: tryme.o getarraypars.tab.o

#getarraypars.tab.o: getarraypars.tab.c getarraypars.tab.h
getarraypars.tab.c getarraypars.tab.h: getarraypars.y

BISON = bison -d

%.tab.c %.tab.h: %.y
	$(BISON) $<

tryme.o: tryme.cpp

local_clean::
	-$(RM) -rf tryme

.PHONY: force


