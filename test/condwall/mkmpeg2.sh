#!/bin/bash
#$Id: mkmpeg2.sh,v 1.1 2005/03/29 21:05:26 mstorti Exp $

yuv_cat() {
    for k in `seq 0 $nstep` 
      do  zyuv=$1/condwall.$k.yuv.gz
      if [ -f $zyuv ] ; then zcat $zyuv ; fi
    done
}

nstep=2000
yuv_cat YUV
