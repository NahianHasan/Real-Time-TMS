function [tri,pp,E_aux,Anor]=genrecipmapks(te2p,p,conductivity,teid,that,rs,ks,omega,scth,th_hair,N,FEMord)
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

%find scalp locations
%define triangle mesh
[tri]=surftri(p',te2p');
if scth~=0;
pcen=(p(:,tri(:,1))+p(:,tri(:,2))+p(:,tri(:,3)))/3;
roicen=sum(p(:,te2p(1,teid)),2)/numel(te2p(1,teid));

err=abs(pcen(1,:)-roicen(1))+abs(pcen(2,:)-roicen(2))+abs(pcen(3,:)-roicen(3));
roicen=pcen(:,err==min(err));
tri=tri(pcen(1,:)-roicen(1)<scth,:);
pcen=pcen(:,pcen(1,:)-roicen(1)<scth);
tri=tri(pcen(1,:)-roicen(1)>-scth,:);
pcen=pcen(:,pcen(1,:)-roicen(1)>-scth);
tri=tri(pcen(2,:)-roicen(2)<scth,:);
pcen=pcen(:,pcen(2,:)-roicen(2)<scth);
tri=tri(pcen(2,:)-roicen(2)>-scth,:);
pcen=pcen(:,pcen(2,:)-roicen(2)>-scth);
tri=tri(pcen(3,:)-roicen(3)<scth,:);
pcen=pcen(:,pcen(3,:)-roicen(3)<scth);
tri=tri(pcen(3,:)-roicen(3)>-scth,:);
pcen=pcen(:,pcen(3,:)-roicen(3)>-scth);
end
[pp,~,tri]=unique(tri(:));
tri=reshape(tri,[numel(tri)/3 3]);
pp=p(:,pp)';
np=numel(pp)/3;
%
 v1=pp(tri(:,2),:)-pp(tri(:,1),:);
 v2=pp(tri(:,3),:)-pp(tri(:,1),:);
 nhat=cross(v1,v2,2);
 for i=1:numel(nhat)/3
 nhat(i,:)=nhat(i,:)/norm(nhat(i,:));
 end
 clear v1 v2 
%extract point normals and local coordinate systems
Anor=zeros([4 4 np]);
v=zeros([3 3]);
nhatp=zeros(size(pp));

for i=1:numel(tri)/3
    v(:,1)=(pp(tri(i,2),:)-pp(tri(i,3),:));
    v(:,1)=v(:,1)/norm(v(:,1));
    v(:,2)=(pp(tri(i,3),:)-pp(tri(i,1),:));
    v(:,2)=v(:,2)/norm(v(:,2));
    v(:,3)=(pp(tri(i,1),:)-pp(tri(i,2),:));
    v(:,3)=v(:,3)/norm(v(:,3));
for j=1:3
nhatp(tri(i,j),:)=nhatp(tri(i,j),:)+...
   acos(-sum(v(:,mod(j,3)+1).*v(:,mod(j+1,3)+1)))...
    *nhat(i,:);
end
end
for i=1:np
nhatp(i,:)=nhatp(i,:)/norm(nhatp(i,:));
end

thatp1=zeros(size(nhatp));
thatp1(:,1)=1-nhatp(:,1).*nhatp(:,1);
thatp1(:,2)=-nhatp(:,1).*nhatp(:,2);
thatp1(:,3)=-nhatp(:,1).*nhatp(:,3);
for i=1:np
thatp1(i,:)=thatp1(i,:)/norm(thatp1(i,:));
end
thatp2=cross(nhatp,thatp1);
pp=pp+nhatp*th_hair;
Anor(1:3,1,:)=reshape(thatp1',[3 1 np]);
Anor(1:3,2,:)=reshape(thatp2',[3 1 np]);
Anor(1:3,3,:)=reshape(nhatp',[3 1 np]);
Anor(1:3,4,:)=reshape(pp',[3 1 np]);
Anor(4,4,:)=1;
clear thatp1 thatp2 nhatp;

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
E_aux=zeros([360 np]);
for j=1:np
for i=1:360
    E_aux(i,j)=-omega*sum(kaux(3,:,i).*(Anor(1:3,3,j)'*Hobs(:,:,j)));
end
end