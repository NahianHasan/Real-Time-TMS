function Neighbors=findneighbors(t2p,p)
edges=[t2p(:,[2,3]);...
       t2p(:,[1,3]);...
       t2p(:,[1,2])];
edges=sort(edges,2);
[edges]=unique(edges,'rows');
edges(end+1:2*end,:)=edges(:,[2,1]);
[~,I] = sort(edges(:,1));
edges=edges(I,:);
Neighbors.list=zeros([numel(edges)/2 1]);
Neighbors.len=zeros([numel(edges)/2 1]);

Neighbors.st=zeros([edges(end,1) 1]);
Neighbors.end=zeros([edges(end,1) 1]);
edid=1;
Neighbors.st(1)=1;
for i=1:numel(edges)/2
    
Neighbors.list(i)=edges(i,2);
Neighbors.len(i)=norm(p(edges(i,1),:)-p(edges(i,2),:));
if edid~=edges(i,1)
Neighbors.end(edid)=i-1;
edid=edid+1;
Neighbors.st(edid)=i;
end

end

if Neighbors.end(end)<Neighbors.st(end)
Neighbors.end(end)=Neighbors.st(end);
end

end