function    Es = spatialEnergy(seg,  width, height, class_num)
seg = reshape(seg, width, height);
[s, t, K] = size(seg);
Xu1=zeros(s,t,K);
Xu1(2:s,2:t,:)=seg(1:s-1,1:t-1,:);

Xu=zeros(s,t,K);
Xu(2:s,:,:)=seg(1:s-1,:,:);

Xur=zeros(s,t,K);
Xur(2:s,1:t-1,:)=seg(1:s-1,2:t,:);

Xr=zeros(s,t,K);
Xr(:,1:t-1,:)=seg(:,2:t,:);

Xdr=zeros(s,t,K);
Xdr(1:s-1,1:t-1,:)=seg(2:s,2:t,:);

Xd=zeros(s,t,K);
Xd(1:s-1,:,:)=seg(2:s,:,:);

Xd1=zeros(s,t,K);
Xd1(1:s-1,2:t,:)=seg(2:s,1:t-1,:);

X1=zeros(s,t,K);
X1(:,2:t,:)=seg(:,1:t-1,:);
% for i=1:class_num
%     XN=sum(cat(3,Xu1==i,Xu==i,Xur==i,Xr==i,Xdr==i,Xd==i,Xd1==i,X1==i),3)/8;
%     Es(:,i) = imstack2vectors(XN);
% end

Es=sum(cat(3,Xu1==seg,Xu==seg,Xur==seg,Xr==seg,Xdr==seg,Xd==seg,Xd1==seg,X1==seg),3)/8; 
Es = imstack2vectors(Es);
Es(:,2) = 1 - Es(:,1);
Es = 1- Es;
end
