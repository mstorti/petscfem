#!/bin/bash
#$Id: mkmpeg2.sh,v 1.9 2003/11/03 03:51:35 mstorti Exp $

nstep=360

yuv_cat() {
    for k in `seq 0 $nstep` 
      do  zyuv=$1/cubcav.$k.yuv.gz
      if [ -f $zyuv ] ; then zcat $zyuv ; fi
    done
}

# yuv_cat YUV-z0
# yuv_cat YUV-z0.25
# yuv_cat YUV-z0.5
# yuv_cat YUV-y0.05
# yuv_cat YUV-x0.95
# yuv_cat YUV-x0.05
# yuv_cat YUV-box

# nstep=600
# yuv_cat YUV-iso
# yuv_cat YUV-iso2
# yuv_cat YUV-iso3

nstep=2000
yuv_cat YUV-streamlines-z05
yuv_cat YUV-streamlines-z01
yuv_cat YUV-limit-streamlines2
yuv_cat YUV-limit-streamlines-below
yuv_cat YUV-core-streamlines
# yuv_cat YUV

#yuv_cat YUV-Re1e6-z0.25
#yuv_cat YUV-Re1e6-z0.01
