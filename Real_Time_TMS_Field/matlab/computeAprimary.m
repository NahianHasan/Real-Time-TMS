function [Aout,curlAout]=computeAprimary(rs,js,nsource,robs,ntarget)
%compute Eobs at observation points (mostly for testing C++ to Fortran interface)
Aout=zeros([3 ntarget]);
curlAout=zeros([3 ntarget]);
iprecEp=10^-4;%FMM accuracy flag lower is faster but less accurate
mex_id_ = 'Aprim(i int, i int, i double[xx], i double[xx], i double[xx], io double[xx], io double[xx], i double)';
[Aout, curlAout] = FEM(mex_id_, nsource, ntarget, rs, js, robs, Aout, curlAout, iprecEp, 3, nsource, 3, nsource, 3, ntarget, 3, ntarget, 3, ntarget);

