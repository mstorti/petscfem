global A B C D bu bp sU sp omega r0 homog_version

omega =0.5;

load -force mat.dat

for j=1:4
  indx=find(mat(:,1)==-1);
endfor

indx=[0;indx];

A=mat((indx(1)+1):(indx(2)-1),:);
B=mat((indx(2)+1):(indx(3)-1),:);
C=mat((indx(3)+1):(indx(4)-1),:);
D=mat((indx(4)+1):(indx(5)-1),:);

clear mat
