function [Hout]=computeHprimary(rs,js,nsource,robs,ntarget)
%compute Eobs at observation points (mostly for testing C++ to Fortran interface)
Hout=zeros([3 ntarget]);
iprecEp=10^-6;%FMM accuracy flag lower is faster but less accurate
mex_id_ = 'Hprim(i int, i int, i double[xx], i double[xx], i double[xx], io double[xx], i double)';
[Hout] = FEM(mex_id_, nsource, ntarget, rs, js, robs, Hout, iprecEp, 3, nsource, 3, nsource, 3, ntarget, 3, ntarget);

