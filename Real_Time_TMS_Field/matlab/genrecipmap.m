function genrecipmap(te2p,p,conductivity,teid,that,rs_e,js_e,rs_h,ks_h,N,FEMord)


%generate auxiliary coils
if numel(rs_e)~=0
Nj=numel(rs_e)/3;
rs2=zeros([3 Nj 360]);
ks2=zeros([3 Nj 360]);
for i=1:360
    phi=i*pi/180;
rs2(3,:,i)=rs_e(3,:);
ks2(:,:,i)=js_e;
rs2(1,:,i)=rs_e(1,:)*cos(phi)+rs_e(2,:)*sin(phi);
rs2(2,:,i)=-rs_e(1,:)*sin(phi)+rs_e(2,:)*cos(phi);
end
[raux_e,jaux_e]=resamplecoil(rs2,ks2,N,360);

clear rs2 ks2;
elseif numel(rs_h)~=0

end
[rv11,jv11]=runcoderecipTensor(te2p,p,conductivity,teid(:)-1,that,FEMord);
