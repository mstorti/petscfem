#__INSERT_LICENSE__
#$Id: Makefile,v 1.2 2001/04/01 03:11:10 mstorti Exp $
SHELL=/bin/sh

.PHONY: save

save: 
	if ( [ -f source.shar,v ] || [ -f RCS/source.shar,v ] ) ; \
		then co -l source.shar ; \
	fi
	if [ ! -f LOG.txt ] ; then touch LOG.txt ; fi
	sshar *.m LOG.txt >source.shar
#	shar *.m Makefile LOG.txt >source.shar
	ci  source.shar
	echo >> LOG.txt
	echo "============================================="  >> LOG.txt
	date >>LOG.txt
	date -u >>LOG.txt
	echo RCS `rlog -h source.shar | grep "^head" ` >> LOG.txt
	echo >> LOG.txt
#	emacsclient LOG.txt


