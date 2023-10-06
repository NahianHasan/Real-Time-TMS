clear all

load ../examples/siminibs0.mat
[t2p,~,nhat]=surftri(p',te2p(:,:)');
area=sqrt(sum(nhat.^2,2))/2;
conductivity=1;
conductivity(reg==1)=0.465;
conductivity(reg==2)=0.01;
conductivity(reg==3)=1.654;
conductivity(reg==4)=0.126;
conductivity(reg==5)=0.276;

cat1=[0.0238;0.04222;-0.001615];
an1=[-0.031841;-0.0160581;-0.0113714];
[~,cat1]=min(...
    (p(1,:)-cat1(1)).^2+ ...
    (p(2,:)-cat1(2)).^2+ ...
    (p(3,:)-cat1(3)).^2 );
[~,an1]=min(...
    (p(1,:)-an1(1)).^2+ ...
    (p(2,:)-an1(2)).^2+ ...
    (p(3,:)-an1(3)).^2 );
rhs=zeros([numel(p)/3 1]);
rhs(cat1)=1;
rhs(an1)=-1;
[x1]=runcodeTES(te2p,p,conductivity,rhs,1);

[scalp,~,t2p]=unique(t2p(:));
t2p=reshape(t2p,[numel(t2p)/3 3]);
scalpp=p(:,scalp)';
%%
trisurf(t2p,scalpp(:,1),scalpp(:,2),scalpp(:,3),(x1(scalp)),'facealpha',.1,'edgealpha',0)
addpath('/media/ljg24/Extreme SSD/neuronmesher')
[pcyl,te2pcyl,~,reg]=loadfem(strcat('/media/ljg24/Extreme SSD/neuronmesher/meshes/cylm',num2str(16),'neuron80.msh'));
rad=1;
A=[rad 0 0 0;0 rad 0 0; 0 0 10 0;0 0 -.0132 1];
pcyl(4,:)=1;
pcyl=A'*pcyl;pcyl=pcyl(1:3,:);
tcyl=surftri(pcyl',te2pcyl');
hold on
trisurf(tcyl,pcyl(1,:)',pcyl(2,:)',pcyl(3,:)','facealpha',1,'edgealpha',0,'facecolor','black')


cat1=[.0377249;.00732004;.000954622];
an1=[-0.0186828;-0.042047;-0.0134186];

[~,cat1]=min(...
    (p(1,:)-cat1(1)).^2+ ...
    (p(2,:)-cat1(2)).^2+ ...
    (p(3,:)-cat1(3)).^2);
[~,an1]=min(...
    (p(1,:)-an1(1)).^2+ ...
    (p(2,:)-an1(2)).^2+ ...
    (p(3,:)-an1(3)).^2 );
rhs=zeros([numel(p)/3 1]);
rhs(cat1)=1;
rhs(an1)=-1;
[x,te2p2,p2,te2te]=runcodeTES(te2p,p,conductivity,rhs,1);
TR = triangulation(te2p',p(1,:)',p(2,:)', p(3,:)');

trisurf(t2p,scalpp(:,1),scalpp(:,2),scalpp(:,3),x(scalp),'edgealpha',0)
save TESheadresult.mat

figure

freq=2000;
beat=5;
amp=5;
duration=.1;
dt=4*10^-5

p=pcyl';te2p=te2pcyl';    
p=p(:,[1 3 2]);
EL=-0.07;
sig=[1 2];
sig=sig(reg);
len=sqrt(sum((p(te2p(:,1),:)-p(te2p(:,2),:)).^2,2));


addpath('D:\neuronmesher\TMS_Efield_Solvers\FEM_MEX_C_codes\');
FEMord=1;
[te2p2,p2]=femgenmesh_c(te2p',p',FEMord);
p=p2';te2p=te2p2'+1;    

[p,te2p,p_mem,t_mem,pb_in,pb_ex,pp_in,pp_ex,pout]=...
    gen_bidomain_mesh(p,te2p,reg,FEMord);
%%%%
   [TI, BC] = pointLocation(TR, p);
xneuron1=sum(BC.*x(TR.ConnectivityList(TI,:)),1);
xneuron2=sum(BC.*x1(TR.ConnectivityList(TI,:)),1);
xneuron1=xneuron1(pout);
xneuron2=xneuron2(pout);
%%build base matrices

laplace=-0.5*femassemble(te2p'-1,p',sig,FEMord);
if FEMord==2
G_mem=grammsurf2ndord(t_mem,p_mem);
else
G_mem=grammsurf(t_mem,p_mem);
end
%%%build t+1 matrix

rl=1;
cm=10^-2;
gl=3;gk=360;gna=1200;
El=-0.0544;Ek=-0.077;Ena=0.05;
Erest=-0.07;
Ntime=round(duration/dt);

soln=1.729810736470345e-01*ones([numel(p)/3,Ntime]);
soln(pp_in,1)= 3.798138360467137e-02;
soln(pp_ex,1)=1.029812168060905e-01;
load  initialcond.mat init;
% soln=init;
rhsF1=zeros([numel(p)/3 1]);
rhsF1(pout)=xneuron1;
rhsF2=zeros([numel(p)/3 1]);
rhsF2(pout)=xneuron2;
tri=surftri(p,te2p);
pcen=(p(tri(:,1),:)+p(tri(:,2),:)+p(tri(:,3),:))/3;
nhat=cross(p(tri(:,2),:)-p(tri(:,1),:),...
           p(tri(:,3),:)-p(tri(:,1),:));

hold on

nmhp=zeros([numel(pb_in) 3]);
nmhp(:,1)=3.176862580749554e-01;
nmhp(:,2)=5.294913493474710e-02; 
nmhp(:,3)=5.963697918631813e-01;


for t=1:Ntime-1
nmh=nmhp;
giter=(gl+gk*nmh(:,1).^4+gna*nmh(:,2).^3.*nmh(:,3))*0.5;
Ein=(gl*El+gk*nmh(:,1).^4*Ek+gna*nmh(:,2).^3.*nmh(:,3)*Ena);


NewTmat(pb_in,pb_in)=-laplace(pb_in,pb_in)+G_mem*(...
    sparse(1:numel(pb_in),1:numel(pb_in),cm/dt+giter));
NewTmat(pb_in,pb_ex)=-laplace(pb_in,pb_ex)+G_mem*(...
    sparse(1:numel(pb_in),1:numel(pb_in),-cm/dt-giter));
NewTmat(pb_ex,pb_in)=-laplace(pb_ex,pb_in)+G_mem*(...
    sparse(1:numel(pb_in),1:numel(pb_in),-cm/dt-giter));
NewTmat(pb_ex,pb_ex)=-laplace(pb_ex,pb_ex)+G_mem*(...
    sparse(1:numel(pb_in),1:numel(pb_in),cm/dt+giter));
%build t matrix
OldTmat(pb_in,pb_in)=laplace(pb_in,pb_in)+G_mem*...
    sparse(1:numel(pb_in),1:numel(pb_in),cm/dt-giter);
OldTmat(pb_in,pb_ex)=laplace(pb_in,pb_ex)+G_mem*...
    sparse(1:numel(pb_in),1:numel(pb_in),-cm/dt+giter);
OldTmat(pb_ex,pb_in)=laplace(pb_ex,pb_in)+G_mem*...
    sparse(1:numel(pb_in),1:numel(pb_in),-cm/dt+giter);
OldTmat(pb_ex,pb_ex)=laplace(pb_ex,pb_ex)+G_mem*...
    sparse(1:numel(pb_in),1:numel(pb_in),cm/dt-giter);

rhs2=rhsF1*sin(2*pi*freq*(dt*(t+1)-.0005))*(dt*t>.0005);%.*(dt*t<1/freq+.0005);%injection current
rhs2=rhsF2*sin(2*pi*(freq+beat)*(dt*(t+1)-.0005))*(dt*t>.0005);%.*(dt*t<1/freq+.0005);%injection current

rhs2(pb_in)=rhs2(pb_in)+G_mem*Ein;
rhs2(pb_ex)=rhs2(pb_ex)-G_mem*Ein;

tic
soln(:,t+1)=NewTmat\(OldTmat*soln(:,t)+rhs2);
toc
if condest(NewTmat)>10^14
soln(:,t+1);
end
nmhp=hhconduct(soln(pb_in,t)-soln(pb_ex,t),nmhp,dt);

 

end
soln(end+1,:)=0;
 trisurf(t_mem,p_mem(:,1),p_mem(:,2),p_mem(:,3),...
     soln(pb_in,t)-soln(pb_ex,t),'facecolor','interp','edgealpha',0)
 caxis([-.08 .05])
 colorbar
%%
col{1}='b'
col{2}='--r'
col{3}='--green'
col{4}='--cyan'
col{5}='--black'
col{6}='--yellow'
for ires=1


tr=triangulation(te2p(:,1:4),p);
y=10^-10:dt:.1;y=y(:);

subplot(4,1,1),
plot(y,...
    amp*(sin(2*pi*freq.*(y-.0005))-sin(2*pi*(freq+beat).*(y-.0005))).*(y>.0005),col{ires})
hold on
ylabel('I (A)')
xlabel('t (s)')
[~,id2]=min(abs(p_mem(:,1))+ ...
    abs(p_mem(:,2)-0.004/4)+ ...
    abs(p_mem(:,3)-10^-6));
subplot(4,1,2),
plot(dt*(1:Ntime),soln(pb_in(id2),:)-soln(pb_ex(id2),:),col{ires})
hold on
ylabel('V (V)')
xlabel('t (s)')

[~,id2]=min(abs(p_mem(:,1))+ ...
    abs(p_mem(:,2)-2*0.004/4)+ ...
    abs(p_mem(:,3)-10^-6));
subplot(4,1,3),
plot(dt*(1:Ntime),soln(pb_in(id2),:)-soln(pb_ex(id2),:),col{ires})
hold on
ylabel('V (V)')
xlabel('t (s)')
[~,id2]=min(abs(p_mem(:,1))+ ...
    abs(p_mem(:,2)-3*0.004/4)+ ...
    abs(p_mem(:,3)+10^-6));
subplot(4,1,4),
plot(dt*(1:Ntime),soln(pb_in(id2),:)-soln(pb_ex(id2),:),col{ires})
hold on
ylabel('V (V)')
xlabel('t (s)')

end


