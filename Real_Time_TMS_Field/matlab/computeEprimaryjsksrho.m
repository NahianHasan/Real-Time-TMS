function [Eout]=computeEprimaryjsksrho(rs,js,ks,rho,nsource,robs,ntarget)
%compute Eobs at observation points (mostly for testing C++ to Fortran interface)
Eout=zeros([3 ntarget]);
iprecEp=10^-4;%FMM accuracy flag lower is faster but less accurate
mex_id_ = 'Eprimrsjsksrho(i int, i int, i double[xx], i double[xx], i double[xx], i double[x], i double[xx], io double[xx], i double)';
[Eout] = FEM(mex_id_, nsource, ntarget, rs, js, ks, rho, robs, Eout, iprecEp, 3, nsource, 3, nsource, 3, nsource, nsource, 3, ntarget, 3, ntarget);


