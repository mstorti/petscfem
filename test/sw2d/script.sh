#!/bin/bash
for n in `seq 0 1 300`
    do echo $n
    convert YUV/plano.$n.tiff YUV/plano.$n.yuv 
    gzip YUV/plano.$n.yuv 
done

## convert  -geometry 800x600 YUV/testsw2d.$n.tiff YUV/testsw2d.$n.tiff 
