function flux=upflux2(U,n);

  global gravity addvisc

  g=gravity;
  t=[n(2) -n(1)];

  NN=size(U,1);
  Uh=(U(2:NN,:)+U(1:NN-1,:))/2;
  h=Uh(:,1);
  sgh=sqrt(gravity*h);
  isgh=1./sgh;
  u=Uh(:,2)./h;
  v=Uh(:,3)./h;

  un=u*n(1)+v*n(2);
  flux=[h.*un h.*un.*u+(g*n(1)/2)*h.^2  h.*un.*v+(g*n(2)/2)*h.^2];
  o=ones(NN-1,1);
  z=zeros(NN-1,1);
  An=[zeros(NN-1,1),-un.*u+g*h*n(1),-un.*v+g*h*n(2), ...
      n(1)*ones(NN-1,1),n(1)*u+un,n(1)*v,n(2)*ones(NN-1,1),n(2)*u,n(2)*v+un];

  d=[un+sgh,un-sgh,un];
  d=abs(d);
  md=max(d')';
  for k=1:3
    d(:,k)=d(:,k)+addvisc*(md-d(:,k));
  end

  %% aA=rowlefs(d,[o,o,z,isgh*n(1),-isgh*n(1),t(1)*o,isgh*n(2),-isgh*n(2),t(2)*o]);  ...
  aA=[d(:,1),d(:,2),z,d(:,1).*isgh*n(1),-d(:,2).*isgh*n(1), ...
d(:,3).*t(1),d(:,1).*isgh*n(2),-d(:,2).*isgh*n(2),d(:,3).*t(2)];

  aA=rowprod([[o,sgh*n(1),sgh*n(2),o,-sgh*n(1), ...
-sgh*n(2)]/2,z,t(1)*o,t(2)*o],aA,3,3,3);

  %%  aA=rowprod(aA,[o -u -v z o z z z o],3,3,3);
  aA(:,1)=aA(:,1)-u.*aA(:,4)-v.*aA(:,7);
  aA(:,2)=aA(:,2)-u.*aA(:,5)-v.*aA(:,8);
  aA(:,3)=aA(:,3)-u.*aA(:,6)-v.*aA(:,9);

  %%  aA=rowprod([o +u +v z o z z z o],aA,3,3,3);
  aA(:,2)=aA(:,2)+u.*aA(:,1);
  aA(:,5)=aA(:,5)+u.*aA(:,4);
  aA(:,8)=aA(:,8)+u.*aA(:,7);

  aA(:,3)=aA(:,3)+v.*aA(:,1);
  aA(:,6)=aA(:,6)+v.*aA(:,4);
  aA(:,9)=aA(:,9)+v.*aA(:,7);

  
  flux=flux-rowprod(aA,(U(2:NN,:)-U(1:NN-1,:))/2,3,1,3);
