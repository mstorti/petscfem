#!/bin/bash
#$Id: mkmpeg2.sh,v 1.4 2003/09/08 10:20:39 mstorti Exp $

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

for k in `seq 0 $nstep` 
do  zcat YUV-y0.05/cubcav.$k.yuv.gz
done

for k in `seq 0 $nstep` 
do  zcat YUV-x0.95/cubcav.$k.yuv.gz
done
