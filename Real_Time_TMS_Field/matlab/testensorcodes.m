clear all
addpath('../../recipcodes/tensorcodes/');

addpath('../../recipcodes/');

%%%%%%% tests reciprocity for each FEM order
mu0=1.25663706*10^-6;
eps0=8.85418782*10^-12;
load ../../recipcodes/lowressphere.mat; 

wm=0.126;csf=1.65;gm=0.276;skull=0.01;skin=.465;air=10^-8;
condu=[wm;gm;csf;skull;skin];
conductivity=zeros(size(reg(:)));
unique(reg(:))'
ct=1;
for cc=unique(reg(:))'%assign conductivities
conductivity(reg==cc)=condu(ct);
ct=ct+1;
end
[~,conductivity]=ndgrid(1:9,conductivity(:));
conductivity([2,3,4,6,7,8],:)=0;

femord=2;
[te2p3,p3,A,a,b,c,d,vol]=FEM_coeffs_neu(p,te2p,conductivity(1,:),femord,1);


%step 1 geometric preprocessing
tic
[te2p2,p2]=femgenmesh_c(te2p,p,femord);
[te2te]=gente2te(te2p);
toc
%% Step 2 assemble FEM matrix

tic
A2=femassembleTensor(te2p2,p2,conductivity,femord);
A3=femassemble(te2p2,p2,conductivity(1,:),femord);

toc
norm(A2(:)-A3(:))/norm(A2(:))


%%
load example1.mat
FEMord=2;
rs=rs';ks=ks';
that1(3,1:181)=that(3,1);
that1(2,1:181)=that(2,1);
that1(1,1:181)=that(1,1);
[~,conductivity]=ndgrid(1:9,conductivity(:));
conductivity([2,3,4,6,7,8],:)=0;

[rv2,jv2]=runcoderecip(te2p,p,conductivity(1,:),teid(:)-1,that1,FEMord);

[rv13,jv13]=runcoderecipTensor(te2p,p,conductivity,teid(:)-1,that1,FEMord);
norm(rv2(:)-rv13(:))/norm(rv13(:))

norm(jv2(:)-jv13(:))/norm(jv13(:))

%%
load example1.mat

[~,conductivity]=ndgrid(1:9,conductivity(:));
conductivity([2,3,4,6,7,8],:)=0;
[E1,x,te2p2,p2,te2te]=runcodeks(te2p,p,conductivity(1,:),rs',ks',rv2,1);%generates Efield at ro1
[E2,x,te2p2,p2,te2te]=runcodekstensor(te2p,p,conductivity,rs',ks',rv2,1);%generates Efield at ro1

p=p';te2p=te2p';
a=p(te2p(:,4),:)-p(te2p(:,1),:);
b=p(te2p(:,3),:)-p(te2p(:,1),:);
c=p(te2p(:,2),:)-p(te2p(:,1),:);
V = abs(1/6 * dot(a', cross(b',c')));
p=p';te2p=te2p';
E_FEM=omega*E1(2,teid)*V(teid)'/sum(V(teid))
imag(H1dotK2_1)
imag(H1dotK2_2)
imag(H1dotK2_3)
