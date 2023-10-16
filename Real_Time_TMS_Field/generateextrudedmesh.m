function [rs,nhat,pp,surf_point_ind,area_tri]=generateextrudedmesh(p,te2p)
    %%%%%%%%generate extruded surface mesh and define inward pointing normals
    [tri,~]=surftri(p',te2p');
    [surf_point_ind,~,tri]=unique(tri(:));
    tri=reshape(tri,[numel(tri)/3 3]);
    pp=p(:,surf_point_ind)';
    v1=pp(tri(:,2),:)-pp(tri(:,1),:);
    v2=pp(tri(:,3),:)-pp(tri(:,1),:);
    nhat=cross(v1,v2,2);

    for i=1:numel(nhat)/3
        nhat(i,:)=nhat(i,:)/norm(nhat(i,:));
    end
    v=zeros([3 3]);
    nhatp=zeros(size(pp));
    for i=1:numel(tri)/3
        v(:,1)=(pp(tri(i,2),:)-pp(tri(i,3),:));
        v(:,1)=v(:,1)/norm(v(:,1));
        v(:,2)=(pp(tri(i,3),:)-pp(tri(i,1),:));
        v(:,2)=v(:,2)/norm(v(:,2));
        v(:,3)=(pp(tri(i,1),:)-pp(tri(i,2),:));
        v(:,3)=v(:,3)/norm(v(:,3));
        for j=1:3
            nhatp(tri(i,j),:)=nhatp(tri(i,j),:)+...
            acos(-sum(v(:,mod(j,3)+1).*v(:,mod(j+1,3)+1)))...
            *nhat(i,:);
        end
    end

    for i=1:numel(pp)/3
        nhatp(i,:)=nhatp(i,:)/norm(nhatp(i,:));
    end
    pp=pp+nhatp*.001;
    %%%%%%%  remove islands  %%%%%%%%%%%%%%%%%
    groups=concomptri(tri,pp);
    group_ids = unique(groups);
    largest_group = group_ids(1);
    largest_count = sum(groups==group_ids(1),'all');
    for ix=2:length(group_ids)
        tx = sum(groups==group_ids(ix),'all');
        if tx > largest_count
            largest_count = tx;
            largest_group = group_ids(ix);
        end
    end
    tri = tri(groups==largest_group,:);
    %%%%%%%%%%%%%%%%%%end surface mesh generation

    v1=pp(tri(:,2),:)-pp(tri(:,1),:);
    v2=pp(tri(:,3),:)-pp(tri(:,1),:);
    rs=(pp(tri(:,1),:)+pp(tri(:,2),:)+pp(tri(:,3),:))'/3;
    nhat=cross(v1,v2,2);
    area_tri = vecnorm(nhat,2,2)/2;
    nhat=nhat'/2;%note nhat is already premultiplied with area
end

function groups=concomptri(t2p,p)
    nt=numel(t2p)/3;
    edges=[t2p(:,[2,3]);
           t2p(:,[1,3]);
           t2p(:,[1,2])];
    node4=[t2p(:,1);t2p(:,2);t2p(:,3)];
    node5=[(1:nt)';(1:nt)';(1:nt)';];
    edges=sort(edges,2);
    [ed,ix,tri2ed]=unique(edges,'rows');
    tri2ed=reshape(tri2ed,[nt 3]);
    ed2tri=zeros([numel(ed)/2 2]);
    ct=1;
    for i=1:numel(tri2ed)
        if ed2tri(tri2ed(i),1)==0
           ed2tri(tri2ed(i),1)=node5(i);
        else
           ed2tri(tri2ed(i),2)=node5(i);
        end
    end

    faces=zeros([nt 1]);

    while nnz(faces==0)~=0
        ids=1:nt;
        ids=ids(faces==0);
        ids=ids(1);
        faces(ids)=ct;
        n=1;
        tids=1;
        while n<=tids
            for edid=1:3
                if ed2tri(tri2ed(ids(n),edid),1)==ids(n)
                    if faces(ed2tri(tri2ed(ids(n),edid),2))==0
                        faces(ed2tri(tri2ed(ids(n),edid),2))=ct;
                        tids=tids+1;
                        ids(tids)=ed2tri(tri2ed(ids(n),edid),2);
                    end
                elseif ed2tri(tri2ed(ids(n),edid),2)==ids(n)
                    if faces(ed2tri(tri2ed(ids(n),edid),1))==0
                        faces(ed2tri(tri2ed(ids(n),edid),1))=ct;
                        tids=tids+1;
                        ids(tids)=ed2tri(tri2ed(ids(n),edid),1);
                    end
                end
            end
            n=n+1;
        end
        nnz(faces==ct);
%         subplot(3,3,ct),
%         hold on
%         trisurf(t2p(faces==ct,:),p(:,1),p(:,2),p(:,3),'edgealpha',0)
        ct=ct+1;
    end
    groups=faces;

end
