#!/bin/bash
launch() { 
    nohup &> nohup-run$1.out make \
		sub_dir=run$1 node_list="$2" run_subdir_test & 
  } 

launch 1 'node14 node15 node17' 
launch 2 'node18 node19'
launch 3 'node20 node21'
launch 4 'node22 node23' 
