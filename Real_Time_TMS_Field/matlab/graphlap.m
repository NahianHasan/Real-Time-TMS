function [Wc,Ac]=graphlap(X,F);

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
cota = @(a,b)cot( acos( dot(normalize(a),normalize(b)) ) );
I = []; J = []; V = []; % indexes to build the sparse matrices
Ia = []; Va = []; % area of vertices
for i=1:3
    % opposite edge e_i indexes
    s = mod(i,3)+1;
    t = mod(i+1,3)+1;
    % adjacent edge
    ctheta = cota(XF(s)-XF(i), XF(t)-XF(i));
    % ctheta = max(ctheta, 1e-2); % avoid degeneracy
    % update the index listing
    I = [I, F(s,:), F(t,:), F(t,:), F(s,:)];
    J = [J, F(t,:), F(s,:), F(t,:), F(s,:)];
    V = [V, ctheta./A, ctheta./A, -ctheta./A, -ctheta./A];
    % update the diagonal with area of face around vertices
    Ia = [Ia, F(i,:)];
    Va = [Va, A];
end
% Aread diagonal matrix
Wc = sparse(I,J,V,n,n)/6;
Ac=sparse(1:n,1:n,ones([n 1]));