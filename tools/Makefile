.PHONY: all
all: odoc.1 odoc.txt

odoc.1: odoc.pl
	pod2man $< > $@

odoc.txt: odoc.pl
	pod2text $< > $@
