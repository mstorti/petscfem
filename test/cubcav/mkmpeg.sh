#!/bin/bash
#$Id: mkmpeg.sh,v 1.1 2003/09/06 22:11:22 mstorti Exp $

mpeg=cubcav.mpeg

if [ -f $mpeg ]
then 
    echo "$mpeg exists. Please delete or rename..."
    exit 1
fi

mkmpeg2.sh | ffmpeg -r 20/1 -s 640x480 -f rawvideo -i - -qscale 5 cubcav.mpeg
