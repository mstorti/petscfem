#!/bin/bash

for k in `seq 0 0`
do  doit.pl $k
    echo before make run
    /usr/bin/make run # > nohup.log
    cat proctable >> nohup.log
    mv nohup.log nohup.log.no-proc$k
done
