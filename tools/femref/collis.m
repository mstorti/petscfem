N=20;				# length of sequence
M=(2:20)';			# different items to be chosen
NN = 1e6;			# Nbr of different sequences 

## Nbr of different possible sequences
Ndseq = stirling(N+M-1)./(stirling(M).*stirling(N-1));

## Estimated nbr of collisions
Ncollis = NN^2 ./Ndseq;

for k = 1:length(M)
  printf("M=%d, Ndseq: %g, Ncollis %g\n",M(k),Ndseq(k),Ncollis(k));
endfor
