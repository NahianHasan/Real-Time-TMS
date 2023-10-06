function [Eout]=FEMinterpolatorks(te2te,te2p,p,rs,js,xval,robs,FEMord)
%compute Etotal at observation points robs
nobs=numel(robs)/3;
np=numel(p)/3;
nsource=numel(rs)/3;
Ep=zeros([3 nobs]);
iprecEp=10^-4;%FMM accuracy flag lower is faster but less accurate
mex_id_ = 'Hprim(i int, i int, i double[xx], i double[xx], i double[xx], io double[xx], i double)';
[Ep] = FEM(mex_id_, nsource, nobs, rs, js, robs, Ep, iprecEp, 3, nsource, 3, nsource, 3, nobs, 3, nobs);
Ep=-Ep;
Eout=zeros([3 nobs]);
if FEMord==1
nte=numel(te2p)/4;
mex_id_ = 'secondaryEfield(i double[xx], i int, i int[xx], i int[xx], i int, i double[xx], i double[x], io double[xx], i double)';
[Eout] = FEM(mex_id_, robs, nobs, te2te, te2p, nte, p, xval, Eout, iprecEp, 3, nobs, 4, nte, 4, nte, 3, np, np, 3, nobs);
elseif FEMord==2
nte=numel(te2p)/10;
mex_id_ = 'secondaryEfield2nd(i double[xx], i int, i int[xx], i int[xx], i int, i double[xx], i double[x], io double[xx], i double)';
[Eout] = FEM(mex_id_, robs, nobs, te2te, te2p, nte, p, xval, Eout, iprecEp, 3, nobs, 4, nte, 10, nte, 3, np, np, 3, nobs);
elseif FEMord==3
nte=numel(te2p)/20;
mex_id_ = 'secondaryEfield3rd(i double[xx], i int, i int[xx], i int[xx], i int, i double[xx], i double[x], io double[xx], i double)';
[Eout] = FEM(mex_id_, robs, nobs, te2te, te2p, nte, p, xval, Eout, iprecEp, 3, nobs, 4, nte, 20, nte, 3, np, np, 3, nobs);
end
Eout=Eout+Ep;


