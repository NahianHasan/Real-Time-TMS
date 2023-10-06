
function phi=findgeodist(X,F,ival,t);

%load mesh
n = size(X,2);
m = size(F,2);



%compute normal and areas
XF = @(i)X(:,F(i,:));
Na = cross( XF(2)-XF(1), XF(3)-XF(1) );
amplitude = @(X)sqrt( sum( X.^2 ) );
A = amplitude(Na)/2;
normalize = @(X)X ./ repmat(amplitude(X), [3 1]);
N = normalize(Na);

%build gradient matrices
I = []; J = []; V = []; % indexes to build the sparse matrices
for i=1:3
    % opposite edge e_i indexes
    s = mod(i,3)+1;
    t = mod(i+1,3)+1;
    % vector N_f^e_i
    wi = cross(XF(t)-XF(s),N); %that
    % update the index listing
    I = [I, 1:m];
    J = [J, F(i,:)];
    V = [V, wi];
end

dA = spdiags(1./(2*A(:)),0,m,m);

GradMat = {};
for k=1:3
    GradMat{k} = dA*sparse(I,J,V(k,:),m,n);
end
Grad = @(u)[GradMat{1}*u, GradMat{2}*u, GradMat{3}*u]';


dAf = spdiags(2*A(:),0,m,m);
DivMat = {GradMat{1}'*dAf, GradMat{2}'*dAf, GradMat{3}'*dAf};


Div = @(q)DivMat{1}*q(1,:)' + DivMat{2}*q(2,:)' + DivMat{3}*q(3,:)';

Delta = DivMat{1}*GradMat{1} + DivMat{2}*GradMat{2} + DivMat{3}*GradMat{3};


%solve for geodesic distances

%solve for geodesic distances
Ac=sparse([F(1,:) F(2,:) F(3,:)],[F(1,:) F(2,:) F(3,:)],[A(:)' A(:)' A(:)']);
Nsoln=numel(ival);
delta = zeros(n,Nsoln);
delta(ival+n*(0:Nsoln-1)) = 1;

u = (Ac+t*Delta)\delta;
for j=1:Nsoln
h = Grad(u(:,j));
h = -normalize(h);
delta(:,j)=Div(h);
end
phi = Delta \ delta;
for j=1:Nsoln
phi(:,j)=phi(:,j)-min(phi(:,j));
end
