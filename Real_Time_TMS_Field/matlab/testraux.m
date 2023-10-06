load('/home/ljg24/Downloads/Simnibtranslation/python/exampleprob.mat','rs2','ks2')
rs=rs2';ks=ks2';
clear rs2 ks2
Nj=360;
%generate auxiliary coils
Ncoil=numel(rs)/3;
rs2=zeros([3 Ncoil Nj]);
ks2=zeros([3 Ncoil Nj]);
for i=1:360
    phi=2*(i/360)*pi;
rs2(3,:,i)=rs(3,:);
ks2(:,:,i)=ks;
rs2(1,:,i)= rs(1,:)*cos(phi)+rs(2,:)*sin(phi);
rs2(2,:,i)=-rs(1,:)*sin(phi)+rs(2,:)*cos(phi);
end
%  cos(phi) sin(phi) 0 0
% -sin(phi) cos(phi) 0 0
%     0        0     1 0
%     0        0     0 1
N=[17 17 2];
[raux,kaux]=resamplecoil(rs2,ks2,N,360);

save('/home/ljg24/Downloads/Simnibtranslation/python/auxdip','raux','kaux','rs','ks','N','Nj');