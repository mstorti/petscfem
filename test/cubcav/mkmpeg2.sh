#!/bin/bash
#$Id: mkmpeg2.sh,v 1.2 2003/09/07 17:16:43 mstorti Exp $

nstep=180

for k in `seq 0 $nstep` 
do  zcat YUV-z0/cubcav.$k.yuv.gz
done

for k in `seq 0 $nstep` 
do  zcat YUV-z0.25/cubcav.$k.yuv.gz
done

for k in `seq 0 $nstep` 
do  zcat YUV-z0.5/cubcav.$k.yuv.gz
done

