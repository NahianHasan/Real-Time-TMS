function [rhs]=femgenrhsjsksrho(te2p,p,reg,rs,js,ks,rho,FEMord)
nte=numel(te2p(1,:));
np=numel(p(1,:));
nsource=numel(rs(1,:));
%computes rhs for 1nd order FEM using rs and js as sources and te2p and p as meshes
%convention is second index is for number of elements. For example, te2p is 4 by nte, and p is 3 by np.
rhs=zeros([np 1]);
%internal accuracy parameter
iprec=10^-4;
if FEMord==1
mex_id_ = 'rhsfuncrsjsksrhos(i int, i int, i int[xx], i double[xx], i double[x], i int, i double[xx], i double[xx], i double[xx], i double[x], io double[x], i double)';
[rhs] = FEM(mex_id_, nte, np, te2p, p, reg, nsource, rs, js, ks, rho, rhs, iprec, 4, nte, 3, np, nte, 3, nsource, 3, nsource, 3, nsource, nsource, np);
elseif FEMord==2
mex_id_ = 'rhsfunc2ndrsjsksrhos(i int, i int, i int[xx], i double[xx], i double[x], i int, i double[xx], i double[xx], i double[xx], i double[x], io double[x], i double)';
[rhs] = FEM(mex_id_, nte, np, te2p, p, reg, nsource, rs, js, ks, rho, rhs, iprec, 10, nte, 3, np, nte, 3, nsource, 3, nsource, 3, nsource, nsource, np);
elseif FEMord==3
mex_id_ = 'rhsfunc3rdrsjsksrhos(i int, i int, i int[xx], i double[xx], i double[x], i int, i double[xx], i double[xx], i double[xx], i double[x], io double[x], i double)';
[rhs] = FEM(mex_id_, nte, np, te2p, p, reg, nsource, rs, js, ks, rho, rhs, iprec, 20, nte, 3, np, nte, 3, nsource, 3, nsource, 3, nsource, nsource, np);
end



