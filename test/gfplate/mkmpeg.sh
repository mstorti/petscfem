#!/bin/bash
#$Id: mkmpeg.sh,v 1.1 2005/03/28 16:42:06 mstorti Exp $

mpeg=tuyere.mpeg

if [ -f $mpeg ] ; then rrm $mpeg ; fi
if [ -f $mpeg ] ; then echo "$mpeg exists" ; exit 1 ; fi

mkmpeg2.sh | ffmpeg -r 16/1 -s 640x480 -f rawvideo -i - -b 4000 $mpeg
