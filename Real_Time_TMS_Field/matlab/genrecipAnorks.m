function [tri,pp,Eaux,Anor]=genrecipAnorks(te2p,p,conductivity,teid,that,rs,ks,omega,Anor,th_hair,N,FEMord)
%generate auxiliary coils
Nj=numel(rs)/3;
rs2=zeros([3 Nj 360]);
ks2=zeros([3 Nj 360]);
for i=1:360
    phi=i*pi/180;
rs2(3,:,i)=rs(3,:);
ks2(:,:,i)=ks;
rs2(1,:,i)= rs(1,:)*cos(phi)+rs(2,:)*sin(phi);
rs2(2,:,i)=-rs(1,:)*sin(phi)+rs(2,:)*cos(phi);
end
%  cos(phi) sin(phi) 0 0
% -sin(phi) cos(phi) 0 0
%     0        0     1 0
%     0        0     0 1
[raux,kaux]=resamplecoil(rs2,ks2,N,360);
clear rs2 ks2;
raux(4,:)=1;
robs=zeros([size(raux) np]);
for i=1:np
robs(:,:,i)=Anor(:,:,i)*raux;
end
robs=reshape(robs(1:3,:,:),[3 prod(N)*np]);

[rv11,jv11]=runcoderecipTensor(te2p,p,conductivity,teid(:)-1,that,FEMord);
tic
[Hobs]=computeHprimary(rv11,jv11,numel(rv11)/3,robs,numel(robs)/3);
toc
Hobs=reshape(Hobs,[3 prod(N),np]);
Eaux=zeros([360 np]);
for j=1:np
for i=1:360
    Eaux(i,j)=-omega*sum(kaux(3,:,i).*(Anor(1:3,3,j)'*Hobs(:,:,j)));
end
end