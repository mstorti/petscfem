#!/bin/bash
#$Id: mkmpeg2.sh,v 1.1 2005/03/28 16:42:06 mstorti Exp $

nstep=360

yuv_cat() {
    for k in `seq 0 $nstep` 
      do  zyuv=$1/tuyere.$k.yuv.gz
      if [ -f $zyuv ] ; then zcat $zyuv ; fi
    done
}


nstep=2000
yuv_cat YUV
