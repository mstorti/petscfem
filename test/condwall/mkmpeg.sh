#!/bin/bash
#$Id: mkmpeg.sh,v 1.2 2005/03/30 01:17:56 mstorti Exp $

mpeg=condwall.mpeg

if [ -f $mpeg ] ; then rrm $mpeg ; fi
if [ -f $mpeg ] ; then echo "$mpeg exists" ; exit 1 ; fi

mkmpeg2.sh | ffmpeg -r 20/1 -s 640x480 -f rawvideo -i - -b 8000 $mpeg
