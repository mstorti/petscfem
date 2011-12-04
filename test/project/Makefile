default: test3d
## default: tests

# BOPT := g_c++

PETSCFEM_DIR := ../..

include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
CASE := bora
CASE_NAME := $(CASE)
PROG := getsurf

tests:  bin test2d test3d

bin:
	$(MAKE) -C ../../tools/project/ project4.bin

reprog:
	$(MAKE) -C ../../tools/project/ -W project4.cpp project4.bin

test2d:
	octave -qH mkdata.m
	$(PETSCFEM_DIR)/tools/project/project4.bin \
		-x xnod1-2D.tmp -i icone1-2D.tmp -s u-2D.tmp \
		-y xnod2-2D.tmp -o ui-2D.tmp -d 2 -e 3 -f 1 
	octave -qH mkcheck.m > verif.tmp

test3d:
	octave -qH mkdata2.m
	$(PETSCFEM_DIR)/tools/project/project4.bin \
		-x xnod1.tmp -i icone1-tet.tmp -s u.tmp \
		-y xnod2.tmp -o ui.tmp -d 3 -e 4 -f 1 
	octave -qH mkcheck2.m > verif3d.tmp

prog:
	$(MAKE) -C $(PETSCFEM_DIR)/tools/project/ project4.bin

link:
	ln -sf $(PETSCFEM_DIR)/tools/project/project4.bin prj

convert:
	$(PETSCFEM_DIR)/tools/hexasplit.bin -i icone1.tmp -o icone1-tet.tmp
