#!/bin/bash
name=movie-doble-rect-1d-ulsar
asf=$name.mpeg

if [ -f $asf ]
then 
    echo "$asf exists. Please delete or rename..."
    exit 1
fi

#mkmpeg2.sh | ffmpeg -r 20/1 -s 1024x768 \
#	-f rawvideo -i - -b 4000 $asf

mkmpeg2.sh | ffmpeg -r 20/1 -s 1200x900 \
	-f rawvideo -i - -qscale 2 $asf
