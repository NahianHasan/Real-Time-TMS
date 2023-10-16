function [Efield,Tm,mem] = real_time_gpu(Q,Ax,X,Y,Z,huygens_surf_points,huygens_surf_points_rot_x,huygens_surf_points_rot_y,huygens_surf_points_rot_z,interp_method,Fields,Fields_temp,SEx,SEy,SEz,SHx,SHy,SHz,Exv,Eyv,Ezv,Hxv,Hyv,Hzv,coeff_g,Efield,Tr,Tri,NModes,memory_calc)
        start_t = tic();
        huygens_surf_points_rot = Tri*huygens_surf_points;
        time_1 = toc(start_t);
        huygens_surf_points_rot_x = huygens_surf_points_rot(1,:);
        huygens_surf_points_rot_y = huygens_surf_points_rot(2,:);
        huygens_surf_points_rot_z = huygens_surf_points_rot(3,:);
        SEx = interpn(X,Y,Z,Exv,huygens_surf_points_rot_x,huygens_surf_points_rot_y,huygens_surf_points_rot_z,interp_method);
        SEy = interpn(X,Y,Z,Eyv,huygens_surf_points_rot_x,huygens_surf_points_rot_y,huygens_surf_points_rot_z,interp_method);
        SEz = interpn(X,Y,Z,Ezv,huygens_surf_points_rot_x,huygens_surf_points_rot_y,huygens_surf_points_rot_z,interp_method);
        SHx = interpn(X,Y,Z,Hxv,huygens_surf_points_rot_x,huygens_surf_points_rot_y,huygens_surf_points_rot_z,interp_method);
        SHy = interpn(X,Y,Z,Hyv,huygens_surf_points_rot_x,huygens_surf_points_rot_y,huygens_surf_points_rot_z,interp_method);
        SHz = interpn(X,Y,Z,Hzv,huygens_surf_points_rot_x,huygens_surf_points_rot_y,huygens_surf_points_rot_z,interp_method);
        Fields_temp = [Tr(1:3,1:3)*[SEx;SEy;SEz];Tr(1:3,1:3)*[SHx;SHy;SHz]./(-4*pi*10^-7)];
        time_2 = toc(start_t);
        for ix=1:length(Fields)
            Fields{ix} = repmat(Fields_temp,1,1,size(Fields{ix},3));
        end
        for ix=1:length(Fields)
            coeff_g{ix} = sum(times(Ax{ix},Fields{ix}),[1,2]);
        end
        coeff_g = squeeze(cat(3,coeff_g{:}));
        time_3 = toc(start_t);
        Efield = Q*coeff_g;
        time_4 = toc(start_t);
        Tm = [time_1,time_2,time_3,time_4];
        mem = [];
        if memory_calc
            mem = memory_calculation(Q,Ax,X,Y,Z,huygens_surf_points,huygens_surf_points_rot_x,huygens_surf_points_rot_y,huygens_surf_points_rot_z,interp_method,Fields,Fields_temp,SEx,SEy,SEz,SHx,SHy,SHz,Exv,Eyv,Ezv,Hxv,Hyv,Hzv,coeff_g,Efield,Tr,Tri,NModes,memory_calc);
        end
end
function [mem] = memory_calculation(varargin)
    mem = 0;
    for ix=1:nargin
        var = gather(varargin{ix});
        mem = mem + ((whos('var').bytes)/1024/1024);
    end
end
