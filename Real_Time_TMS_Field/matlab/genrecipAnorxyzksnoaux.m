function [Eaux]=genrecipAnorxyzksnoaux(te2p,p,conductivity,teid,rs,ks,omega,Anor,N,FEMord)
%generate auxiliary coils
np=numel(Anor(1,1,:));
rs(4,:)=1;
Ns=numel(rs)/4;
robs=zeros([size(rs) np]);
for i=1:np
robs(:,:,i)=Anor(:,:,i)*rs;
end
robs=reshape(robs(1:3,:,:),[3 Ns*np]);
Eaux=zeros([np,3]);
for dir=1:3

    that=zeros([3 numel(teid)]);
    that(dir,:)=1;
[rv11,jv11]=runcoderecipTensor(te2p,p,conductivity,teid(:)-1,that,FEMord);
tic
[Hobs]=computeHprimary(rv11,jv11,numel(rv11)/3,robs,numel(robs)/3);
toc
Hobs=reshape(Hobs,[3 Ns np]);
for j=1:np
    Eaux(j,dir)=-omega*sum(ks(3,:).*(Anor(1:3,3,j)'*Hobs(:,:,j)));
end

end

end

