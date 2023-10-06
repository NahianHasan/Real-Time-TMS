function [Efield_cpu,Tm] = real_time_cpu(Q,Ax,X,Y,Z,huygens_surf_points,interp_method,Exv,Eyv,Ezv,Hxv,Hyv,Hzv,Transformation,NModes)
    start_t = tic();
    huygens_surf_points_rot = Transformation(:,:)\huygens_surf_points;
    time_1 = toc(start_t);
    SEx = interpn(X,Y,Z,Exv,huygens_surf_points_rot(1,:),huygens_surf_points_rot(2,:),huygens_surf_points_rot(3,:),interp_method);
    SEy = interpn(X,Y,Z,Eyv,huygens_surf_points_rot(1,:),huygens_surf_points_rot(2,:),huygens_surf_points_rot(3,:),interp_method);
    SEz = interpn(X,Y,Z,Ezv,huygens_surf_points_rot(1,:),huygens_surf_points_rot(2,:),huygens_surf_points_rot(3,:),interp_method);
    SHx = interpn(X,Y,Z,Hxv,huygens_surf_points_rot(1,:),huygens_surf_points_rot(2,:),huygens_surf_points_rot(3,:),interp_method);
    SHy = interpn(X,Y,Z,Hyv,huygens_surf_points_rot(1,:),huygens_surf_points_rot(2,:),huygens_surf_points_rot(3,:),interp_method);
    SHz = interpn(X,Y,Z,Hzv,huygens_surf_points_rot(1,:),huygens_surf_points_rot(2,:),huygens_surf_points_rot(3,:),interp_method);
    Fields = [Transformation(1:3,1:3)*[SEx;SEy;SEz];Transformation(1:3,1:3)*[SHx;SHy;SHz]./(-4*pi*10^-7)];
    Fields(isnan(Fields)) = 0;
    Fields = reshape(Fields,[],1);
    time_2 = toc(start_t);
    coeff_cpu=zeros([NModes 1]);
    for i=1:NModes
        curr = squeeze(Ax(1:6,:,i));
        coeff_cpu(i)=reshape(curr,1,[])*Fields;
    end
    time_3 = toc(start_t);
    Efield_cpu=Q*coeff_cpu;
    time_4 = toc(start_t);                
    Tm = [time_1,time_2,time_3,time_4];
end