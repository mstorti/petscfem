# $Id: Makefile,v 1.1 2003/02/17 14:53:21 mstorti Exp $ 

.PHONY: $(plano_list)
plano_list = plano_cons plano_cons_auto plano_local plano_auto plano_local_nm \
		plano_local_weak plano_fm2_weak 

#w Test for shallow water
tests: adv $(plano_list)

plano.depl: force

plano_cons:       EPERLFLAGS = "-d consistent_supg=1"
plano_cons_auto:  EPERLFLAGS = "-d consistent_supg_auto=1"
plano_local:      EPERLFLAGS = "-d local_time_step=1"
plano_auto:       EPERLFLAGS = "-d auto_time_step=1"
plano_local_nm:   EPERLFLAGS = "-d local_time_step_nm=1"
plano_local_weak: EPERLFLAGS = "-d local_time_step_weak=1"
plano_fm2_weak:   EPERLFLAGS = "-d fm2_weak=1"

plano_%: 
	$(MAKE) EPERLFLAGS=$(EPERLFLAGS) plano.depl
	-$(MPIRUN) -np 1 $(ADV) -case plano.depl &> $@.sal

local_clean::
	-rm -rf *.sal

SRCS = 
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

#s
