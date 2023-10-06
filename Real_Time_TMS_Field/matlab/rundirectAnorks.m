function [Eaux]=rundirectAnorks(te2p,p,conductivity,teid,rs,ks,omega,Anor,FEMord)
%generate auxiliary coils
rs(4,:)=1;
rs=Anor(1:3,:)*rs;
ks=omega*Anor(1:3,1:3)*ks;
ro=(p(:,te2p(1,teid))+p(:,te2p(2,teid))+p(:,te2p(3,teid))+p(:,te2p(4,teid)))/4;
[Eaux]=runcodekstensor(te2p,p,conductivity,rs,ks,ro,FEMord);