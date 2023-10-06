load  ../examples/example1.mat
rs=rs';ks=ks';
that1(3,1:181)=that(3,1);
that1(2,1:181)=that(2,1);
that1(1,1:181)=that(1,1);
[rv11,jv11]=runcoderecip(te2p,p,conductivity,teid(:)-1,that1,1);
norm(jv11(:)-jv1(:))/norm(jv1(:))
norm(rv11(:)-rv2(:))/norm(rv2(:))
%[rv12,jv12]=runcoderecip(te2p,p,conductivity,teid(:)-1,that1,2);
%norm(jv12(:)-jv2(:))/norm(jv2(:))
%norm(rv12(:)-rv2(:))/norm(rv2(:))
%[rv13,jv13]=runcoderecip(te2p,p,conductivity,teid(:)-1,that1,3);
%norm(jv13(:)-jv3(:))/norm(jv3(:))
%norm(rv13(:)-rv2(:))/norm(rv2(:))
%%
tic
[Hobs]=computeHprimary(rv11(:,1:4000),jv11(:,1:4000),4000,cat(2,rs,rs),2*696);
toc
%%
H1K2_1=omega*sum(sum(Hobs.*ks))
H1dotK2_1
[Hobs]=computeHprimary(rv12,jv12,numel(rv12)/3,rs,696);
H1K2_2=omega*sum(sum(Hobs.*ks))
H1dotK2_2
[Hobs]=computeHprimary(rv13,jv13,numel(rv13)/3,rs,696);
H1K2_3=omega*sum(sum(Hobs.*ks))
H1dotK2_3
%
N=[7 7 2];
Nj=numel(rs)/3;
rs2=zeros([3 Nj 360]);
ks2=zeros([3 Nj 360]);

for i=1:360
    phi=i*pi/180;
rs2(3,:,i)=rs(3,:);
ks2(:,:,i)=ks;
rs2(1,:,i)=rs(1,:)*cos(phi)+rs(2,:)*sin(phi);
rs2(2,:,i)=-rs(1,:)*sin(phi)+rs(2,:)*cos(phi);
end
tic
[Hobs]=computeHprimary(rv13,jv13,numel(rv13)/3,reshape(rs2,[3,360*Nj]),360*Nj);
Hobs=reshape(Hobs,[3,Nj,360]);
Recip1=sum(sum(Hobs.*ks,1),2);
toc

[raux,kaux]=resamplecoil(rs2,ks2,N,360);
tic
[Hobs]=computeHprimary(rv13,jv13,numel(rv13)/3,raux,prod(N));
for i=1:360
Recip2(i)=sum(sum(Hobs(:,:).*kaux(:,:,i)));
end
toc
prod(N)
plot(Recip1(:))
hold on
plot(Recip2(:))
