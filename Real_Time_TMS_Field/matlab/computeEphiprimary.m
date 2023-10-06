function [Ephi]=computeEphiprimary(rs,rho,nsource,robs,ntarget)
%compute Eobs at observation points (mostly for testing C++ to Fortran interface)
Ephi=zeros([3 ntarget]);
iprecEp=10^-4;%FMM accuracy flag lower is faster but less accurate
mex_id_ = 'Ephiprim(i int, i int, i double[xx], i double[x], i double[xx], io double[xx], i double)';
[Ephi] = FEM(mex_id_, nsource, ntarget, rs, rho, robs, Ephi, iprecEp, 3, nsource, nsource, 3, ntarget, 3, ntarget);



