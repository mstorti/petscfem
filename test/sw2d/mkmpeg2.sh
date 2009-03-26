#!/bin/bash

yuv_cat() {
    for k in `seq 0 1 300`
      do  zyuv=YUV/plano.$k.yuv.gz
      if [ -f $zyuv ] ; then zcat $zyuv ; fi
    done

#     for k in `seq 262 1065`
#       do  zyuv=frames/furn_presion.$k.yuv.gz
#       if [ -f $zyuv ] ; then zcat $zyuv ; fi
#     done
}

yuv_cat VIDEO
