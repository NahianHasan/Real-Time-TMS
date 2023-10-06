function [Eaux]=genrecipAnorxyzks(te2p,p,conductivity,teid,rs,ks,omega,Anor,N,FEMord)
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
[raux,kaux]=resamplecoil(rs2,ks2,N,360);
clear rs2 ks2;

np=size(Anor,3);
raux(4,:)=1;
robs=zeros([size(raux) np]);
for i=1:np
robs(:,:,i)=Anor(:,:,i)*raux;
end
robs=reshape(robs(1:3,:,:),[3 prod(N)*np]);
Eaux=zeros([360,np,3]);
for dir=1:3

    that=zeros([3 numel(teid)]);
    that(dir,:)=1;
[rv11,jv11]=runcoderecipTensor(te2p,p,conductivity,teid(:)-1,that,FEMord);
save(strcat('volumecurrs',num2str(dir),'.mat'),'rv11','jv11');
tic
[Hobs]=computeHprimary(rv11,jv11,numel(rv11)/3,robs,numel(robs)/3);
toc
Hobs=reshape(Hobs,[3 prod(N),np]);
save(strcat('bookmark',num2str(dir),'.mat'),'jv11','rv11','Hobs');
for j=1:np
for i=1:360
    Eaux(i,j,dir)=-omega*sum(kaux(3,:,i).*(Anor(1:3,3,j)'*Hobs(:,:,j)));
end
end

end

end

