#!/bin/bash
#$Id: mkmpeg.sh,v 1.2 2003/11/03 03:51:35 mstorti Exp $

mpeg=cubcav.avi

if [ -f $mpeg ] ; then rrm $mpeg ; fi
if [ -f $mpeg ] ; then echo "$mpeg exists" ; exit 1 ; fi

mkmpeg2.sh | ffmpeg -r 20/1 -s 640x480 -f rawvideo -i - -b 2000 $mpeg
