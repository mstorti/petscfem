PETSCFEM_DIR = ..

SRCS = tryme.cpp getarraypars.tab.c

#getarraypars.tab.o: override BOPT = g
include $(PETSCFEM_DIR)/Makefile.base

tryme: LINK.o = g++
tryme: tryme.o getarraypars.tab.o

BISON = bison -d

local_depend: getarraypars.tab.c getarraypars.tab.h

.INTERMEDIATE: *.tab.*
%.tab.c %.tab.h: %.y
	$(BISON) $<

tryme.o: tryme.cpp

%.o: %.c
	$(GNU_C_COMPILER) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

local_clean::
	-$(RM) -rf tryme *.tab.*

