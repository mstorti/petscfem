rm -f output.graph.tmp
f () { 
    echo ================================= >> output.graph.tmp 
    echo N $1 M $2 npart $3 >> output.graph.tmp 
    tryme.bin $1 $2 $3 1 >> output.graph.tmp 
} 
f 1 60 6 1 
f 1 60 7 1 
f 1 60 7 1 
f 1 60 8 1 
f 1 60 9 1 
f 1 60 10 1 
f 1 60 10 1 
f 1 60 10 1 
f 10 60 10 1 
f 10 10 10 1 
f 20 20 10 1
