#!/bin/bash
#$Id: mkmpeg2.sh,v 1.1 2003/09/06 22:11:22 mstorti Exp $

for k in `seq 0 100` 
do  zcat YUV-z0.5/cubcav.$k.yuv.gz
done

for k in `seq 0 100` 
do  zcat YUV-z0.25/cubcav.$k.yuv.gz
done
