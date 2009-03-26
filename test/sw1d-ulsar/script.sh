#!/bin/bash
for n in `seq 0 1 170`
    do echo $n
    convert YUV/doble-rect.$n.tiff YUV/doble-rect.$n.yuv 
    gzip YUV/doble-rect.$n.yuv 
done

## convert  -geometry 800x600 YUV/testsw2d.$n.tiff YUV/testsw2d.$n.tiff 
