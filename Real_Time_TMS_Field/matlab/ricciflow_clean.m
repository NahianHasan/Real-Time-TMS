function [lv,uval,Kmetric,AreaC]=ricciflow_clean(t2p,p);

nt=numel(t2p)/3;
np=numel(p)/3;
edges=[t2p(:,[2,3]);
       t2p(:,[3,1]);
       t2p(:,[1,2])];
edges=sort(edges,2);
[edarr,ix,edid]=unique(edges,'rows','stable');   
vec=hist(edid,1:max(edid));
bdr=edarr(vec==1,:);
bdr=unique(bdr(:));
defect=2*pi*ones([np 1]);
defect(bdr)=pi;
interior=1:np;
interior=interior(defect==2*pi);

%%%%%compute mesh edge lengths
lv(:,1)=sqrt(sum((p(t2p(:,2),:)-p(t2p(:,3),:)).^2,2));
lv(:,2)=sqrt(sum((p(t2p(:,1),:)-p(t2p(:,3),:)).^2,2));
lv(:,3)=sqrt(sum((p(t2p(:,1),:)-p(t2p(:,2),:)).^2,2));

gamma_f=repmat(sum(lv,2),[1 3])-2*lv;
gamma=10^20*ones([np 1]);
for i=1:3*nt
gamma(t2p(i))=min(gamma(t2p(i)),gamma_f(i));
end
clear ltot;
uval=log(gamma);
sql=lv.^2;
gamma=gamma/2;
gamma_f=gamma(t2p);
sqgamma=gamma_f.^2;
%%%%%%%compute invariant indicator function
phi=(sql-sqgamma(:,[2;3;1])-sqgamma(:,[3;1;2]))*.5;
phi=phi./(gamma_f(:,[2;3;1]).*gamma_f(:,[3;1;2]));
clear sqgamma gamma_f
%%%%compute angle defects
theta=(repmat(sum(sql,2),[1 3])-2*sql)*0.5;
theta=acos(theta./(lv(:,[2;3;1]).*lv(:,[3;1;2])));
clear sql;

Area(:,1)=abs(lv(:,2).*lv(:,3).*sin(theta(:,1)))*.5;
Kmetric=zeros([np 1]);
AreaC=zeros([np 2]);
for i=1:3*nt
Kmetric(t2p(i))=Kmetric(t2p(i))-theta(i);
end
for i=1:nt
AreaC(t2p(i,1),1)=AreaC(t2p(i,1),1)+Area(i,1);
AreaC(t2p(i,2),1)=AreaC(t2p(i,2),1)+Area(i,1);
AreaC(t2p(i,3),1)=AreaC(t2p(i,3),1)+Area(i,1);
end
Kmetric=(defect+Kmetric);
GAUSSBONNETNUMBER=sum(Kmetric)

Ktarget=zeros([np 1]);
Ktarget(bdr)=2*pi/numel(bdr);
%%%%%%%%%%%%start algorthm
clear sql
tau=10^-5;
eps=.01;
counter=0;
freebdr=1;
error=1000;
while error>tau && counter<=100
counter=counter+1;
[lv,gamma,Kmetric,theta]=...
    updatetop(uval,phi,t2p,defect,nt,np);
tic
Hess=computehessian(lv,gamma,t2p,phi);
toc
tic
Gamma_f=gamma(t2p); 
if freebdr
Hess(bdr,:) = [];
Hess(:,bdr) = [];
dU=Hess\(Ktarget(interior)-Kmetric(interior));
%dU=eps*(Ktarget-Kmetric);
    uval(interior)=uval(interior)+dU;
 er(counter)=max(abs(Ktarget(interior)-Kmetric(interior)));

else
dU=Hess\(Ktarget(:)-Kmetric(:));
%dU=eps*(Ktarget-Kmetric);
    uval(:)=uval(:)+dU;
er(counter)=max(abs(Ktarget(:)-Kmetric(:)));
end
s=sum(uval(:));
uval=uval-s/np;
error=er(counter);
end
plot(er)
Area(:,2)=abs(lv(:,2).*lv(:,3).*sin(theta(:,1)))*.5;
for i=1:nt
AreaC(t2p(i,1),2)=AreaC(t2p(i,1),2)+Area(i,2);
AreaC(t2p(i,2),2)=AreaC(t2p(i,2),2)+Area(i,2);
AreaC(t2p(i,3),2)=AreaC(t2p(i,3),2)+Area(i,2);
end

end

function Hess=computehessian(lv,gamma,t2p,I);
gamma=gamma(t2p);
nt=numel(t2p)/3;
col=zeros([nt 9]);
row=zeros([nt 9]);
val=zeros([nt 9]);
i=[1;2;3];
j=[2;3;1];
k=[3;1;2];
jk=[1;2;3];
ki=[2;3;1];
ij=[3;1;2];

A=2*lv(:,ij).*lv(:,ki);
B=lv(:,ij).^2+lv(:,ki).^2-lv(:,jk).^2;
C=2*(gamma(:,i)+gamma(:,j).*I(:,ij)).*lv(:,ki)./lv(:,ij)+...
  2*(gamma(:,i)+gamma(:,k).*I(:,ki)).*lv(:,ij)./lv(:,ki);
D=2*(2*gamma(:,i)+gamma(:,j).*I(:,ij)+...
                 +gamma(:,k).*I(:,ki));
denominator=A.*sqrt(A.^2-B.^2);
A=A./denominator;
B=B./denominator;
clear denominator;
col(:,1:3)=t2p(:,i);
row(:,1:3)=t2p(:,i);
val(:,1:3)=A.*D-B.*C;
val(:,1:3)=val(:,1:3).*gamma(:,i);
clear C D;
E=2*(gamma(:,j)+gamma(:,i).*I(:,ij)).*lv(:,ki)./lv(:,ij);
F=2*((gamma(:,i).*I(:,ij)-gamma(:,k).*I(:,jk)));
col(:,4:6)=t2p(:,i);
row(:,4:6)=t2p(:,j);
val(:,4:6)=A.*F-B.*E;
val(:,4:6)=val(:,4:6).*gamma(:,j);

k=[2;3;1];
j=[3;1;2];
ki=[3;1;2];
ij=[2;3;1];

E=2*(gamma(:,j)+gamma(:,i).*I(:,ij)).*lv(:,ki)./lv(:,ij);
F=2*((gamma(:,i).*I(:,ij)-gamma(:,k).*I(:,jk)));
col(:,7:9)=t2p(:,i);
row(:,7:9)=t2p(:,j);
val(:,7:9)=A.*F-B.*E;
val(:,7:9)=val(:,7:9).*gamma(:,j);

Hess=sparse(col(:),row(:),val(:));



end

function [l_new,gamma_new,K_new,theta]=...
    updatetop(uval,phi,t2p,defect,nt,np);
K_new=zeros([np 1]);
gamma_new=exp(uval);
gamma_f=gamma_new(t2p);
sqg=gamma_f.^2;
l_new=sqg(:,[2;3;1])+sqg(:,[3;1;2])+...
    2*gamma_f(:,[2;3;1]).*gamma_f(:,[3;1;2]).*phi;

theta=(repmat(sum(l_new,2),[1 3])-2*l_new)*0.5;
l_new=sqrt(l_new);
theta=acos(theta./(l_new(:,[2;3;1]).*l_new(:,[3;1;2])));

for i=1:3*nt
K_new(t2p(i))=K_new(t2p(i))-theta(i);
end
K_new=(defect(:)+K_new(:));
end