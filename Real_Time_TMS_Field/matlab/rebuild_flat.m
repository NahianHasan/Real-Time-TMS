function [pnew,t2pnew]=rebuild_flat(t2p,p,l_new);
nt=numel(t2p)/3;
np=numel(p)/3;
edges=[t2p(:,[2,3]);
       t2p(:,[3,1]);
       t2p(:,[1,2])];
edges=sort(edges,2);
[edarr,ix,edid]=unique(edges,'rows','stable');   

t2e=reshape(edid,[nt 3]);
e2t=zeros([numel(edarr)/2 2]);
for j=1:3
for i=1:nt
    if e2t(t2e(i,j),1)==0
        e2t(t2e(i,j),1)=i;
    else
        e2t(t2e(i,j),2)=i;
    end        
end
end
neighbors=zeros([nt 3]);
for j=1:3
for i=1:nt
    if e2t(t2e(i,j),1)==i
    neighbors(i,j)=e2t(t2e(i,j),2);
    else
    neighbors(i,j)=e2t(t2e(i,j),1);
    end        
end
end
clear e2t t2e
sql=l_new.^2;
ssql=sum(sql,2);
theta_new=zeros([nt 3]);
for i=1:nt
theta_new(i,1)=acos((ssql(i)-2*sql(i,1))/...
    (2*l_new(i,2)*l_new(i,3)));
theta_new(i,2)=acos((ssql(i)-2*sql(i,2))/...
    (2*l_new(i,1)*l_new(i,3)));
theta_new(i,3)=acos((ssql(i)-2*sql(i,3))/...
    (2*l_new(i,1)*l_new(i,2)));
end

queue=zeros([nt 1]);
queuef=zeros([nt 1]);
include=zeros([nt 1]);
embedded=zeros([np 1]);
embeddedt=zeros([nt 1]);
pnew=zeros([np 3]);
pnew(t2p(1,1),:)=[0,0,0];
pnew(t2p(1,2),:)=[l_new(1,3),0,0];
pnew(t2p(1,3),:)=[l_new(1,2)*cos(theta_new(1,1)),l_new(1,2)*sin(theta_new(1,1)),0];
embedded(t2p(1,:))=1;
embeddedt(1)=1;
ct=0;
for i=1:3
    if neighbors(1,i)~=0
    if embeddedt(neighbors(1,i))==0 
queue(ct+1)=neighbors(1,i);
queuef(ct+1)=t2p(1,i);
ct=ct+1;
    end
    end
end
include(1)=1;
while ct~=0
tr=queue(1);
queue(1)=[];
fr=queuef(1);
queuef(1)=[];
if tr==0 || embeddedt(tr)==1
    ct=ct-1;
    continue;
end

    embeddedt(tr)=1;
if embedded(t2p(tr,1))~=1
    embedded(t2p(tr,1))=1;
    include(tr)=1;
pnew(t2p(tr,1),:)=findpoint(pnew(t2p(tr,2),:),pnew(t2p(tr,3),:),...
    l_new(tr,3),l_new(tr,2),pnew(fr,:));

 elseif embedded(t2p(tr,2))~=1
    embedded(t2p(tr,2))=1;
    include(tr)=1;
pnew(t2p(tr,2),:)=findpoint(pnew(t2p(tr,3),:),pnew(t2p(tr,1),:),...
    l_new(tr,1),l_new(tr,3),pnew(fr,:));

 elseif embedded(t2p(tr,3))~=1
    embedded(t2p(tr,3))=1; 
    include(tr)=1;
pnew(t2p(tr,3),:)=findpoint(pnew(t2p(tr,1),:),pnew(t2p(tr,2),:),...
    l_new(tr,2),l_new(tr,1),pnew(fr,:));
end

ct=ct-1;
for i=1:3
    if neighbors(tr,i)~=0
    if embeddedt(neighbors(tr,i))==0 
queue(ct+1)=neighbors(tr,i);
queuef(ct+1)=t2p(tr,i);
ct=ct+1;
    end
    end
end


    
end
t2pnew=t2p(include==1,:);
end





function p3=findpoint(p1,p2,l1,l2,pn)

    [xout,yout] = circcirc(p1(1),p1(2),l1,p2(1),p2(2),l2);
a1=sum(([xout(1) yout(1)]-pn(1:2)).^2);
a2=sum(([xout(2) yout(2)]-pn(1:2)).^2);
if a1>a2
p3=[xout(1) yout(1) 0];
else
p3=[xout(2) yout(2) 0];
end

end
