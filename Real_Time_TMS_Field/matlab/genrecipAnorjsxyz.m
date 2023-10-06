function [Eaux,Anor]=genrecipAnorjsxyz(te2p,p,conductivity,teid,rs,js,omega,Anor,th_hair,N,FEMord)
%generate auxiliary coils
np=numel(Anor)/16;
Nj=numel(rs)/3;
rs2=zeros([3 Nj 360]);
js2=zeros([3 Nj 360]);
for i=1:360
    phi=i*pi/180;
rs2(3,:,i)=rs(3,:);
js2(:,:,i)=js;
rs2(1,:,i)= rs(1,:)*cos(phi)+rs(2,:)*sin(phi);
rs2(2,:,i)=-rs(1,:)*sin(phi)+rs(2,:)*cos(phi);
end
%  cos(phi) sin(phi) 0 0
% -sin(phi) cos(phi) 0 0
%     0        0     1 0
%     0        0     0 1
[raux,jaux]=resamplecoiljs(rs2,js2,N,360);
clear rs2 ks2;
raux(4,:)=1;
robs=zeros([size(raux) np]);
for i=1:np
robs(:,:,i)=Anor(:,:,i)*raux;
end
robs=reshape(robs(1:3,:,:),[3 prod(N)*np]);
Eaux=zeros([360 np 3]);
for dir=1:3
   that=zeros([3 numel(teid)]); 
   that(:,dir)=1;
[rv11,jv11]=runcoderecipTensor(te2p,p,conductivity,teid(:)-1,that,FEMord);
tic
[Eobs]=computeEprimary(rv11,jv11,numel(rv11)/3,robs,numel(robs)/3);
toc
Eobs=reshape(Eobs,[3 prod(N),np]);
for j=1:np
    Etemp=(Anor(1:3,1:3,j)'*Eobs(:,:,j));
for i=1:360
    Eaux(i,j,dir)=omega*sum(sum(jaux(:,:,i).*Etemp));
end
end
end