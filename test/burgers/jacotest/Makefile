#__INSERT_LICENSE__

.PHONY: force all run tests data

default: tests

tests: test_jaco_1elem

PETSCFEM_DIR := ../../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 

PROG := $(ADVDIF)
LOC_PROG := advdif

test_jaco_1elem: 
	rm -f mat.output *.res
	make-depl burgers.epl burgers.depl -d check_jac=1
	-$(MPIEXEC) -n 1 $(ADVDIF) -case burgers.depl \
		</dev/null &> advdif.log
	$(OCTAVE) proc.m

local_clean::
	rm -f mat.output STEPS*

#w Makes local link for running with debugger
link:
	$(MAKE) -C $(PETSCFEM_DIR) $(LOC_PROG)
	ln -f $(PROG) $(LOC_PROG)

prog:
	$(MAKE) -C $(PETSCFEM_DIR) $(LOC_PROG)
#s
