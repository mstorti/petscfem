#!/bin/bash
#$Id: mkmpeg2.sh,v 1.5 2003/09/09 11:45:11 mstorti Exp $

nstep=360

yuv_cat() {
    for k in `seq 0 $nstep` 
      do  zyuv=$1/cubcav.$k.yuv.gz
      if [ -f $zyuv ] ; then zcat $zyuv ; fi
    done
}

yuv_cat YUV-z0
yuv_cat YUV-z0.25
yuv_cat YUV-z0.5
yuv_cat YUV-y0.05
yuv_cat YUV-x0.95
yuv_cat YUV-x0.05
yuv_cat YUV-box

nstep=600
yuv_cat YUV-iso
yuv_cat YUV-iso2
yuv_cat YUV
