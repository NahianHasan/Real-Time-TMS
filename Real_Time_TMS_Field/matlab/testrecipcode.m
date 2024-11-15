clear all
addpath C:\Users\Admin\Desktop\FEM_modes_new_source_Files\matlab
addpath C:\Users\Admin\Downloads

[vv,ff] = icosphere(4);
vv(:,1)=1.5*vv(:,1);
v1=vv(ff(:,2),:)-vv(ff(:,1),:);
v2=vv(ff(:,3),:)-vv(ff(:,1),:);
rr=(vv(ff(:,1),:)+vv(ff(:,2),:)+vv(ff(:,3),:))/3;
nhat=cross(v1,v2,2);
nhat=nhat/2;
Nsour=100;
rs(:,1)=.032*cos(0:2*pi/Nsour:2*pi);
rs(:,2)=.032*sin(0:2*pi/Nsour:2*pi);
rs(:,3)=0;

js=rs(2:end,:)-rs(1:end-1,:);
rs=(rs(2:end,:)+rs(1:end-1,:))/2;

[Eout]=computeEprimary(rs',js',numel(rs)/3,rr',numel(rr)/3);
[Hout]=computeHprimary(rs',js',numel(rs)/3,rr',numel(rr)/3);
Hout=Hout/(4*pi*10^-7);%curl of A

%%
Ks=-cross(nhat,Eout');
Js=cross(nhat,Hout');
rho=sum(nhat.*(Eout'),2);

robservation=1.1*vv;
[Eout]=computeEprimary(rs',js',numel(rs)/3,robservation',numel(robservation)/3);


[Aprim]=computeEprimary(rr',Js',numel(Js)/3,robservation',numel(robservation)/3);
[Ephi]=computeEphiprimary(rr',rho,numel(Js)/3,robservation',numel(robservation)/3);

%[Ephi]=computeEphiprimary([0 0 0]',1,1,robservation',numel(robservation)/3);


[Fk]=computeHprimary(rr',Ks',numel(Js)/3,robservation',numel(robservation)/3);
Fk=-Fk/(4*pi*10^-7);
Eout2=Aprim+Fk+Ephi;

subplot(3,1,1),
trisurf(ff,vv(:,1),vv(:,2),vv(:,3),Eout(1,:)');
subplot(3,1,2),
trisurf(ff,vv(:,1),vv(:,2),vv(:,3),Eout2(1,:)');
subplot(3,1,3),
plot(Eout(:))
hold on
plot(Eout2(:),'--')


