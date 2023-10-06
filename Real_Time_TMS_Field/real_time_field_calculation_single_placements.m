function [Efield] = real_time_field_calculation_single_placements(Q,Ax,params,Transformation)
    multiplying_factor = 6.6E7;
    if numel(params{31})>0
        Tr = gpuArray(single(Transformation));
        Tri = gpuArray(single(inv(Tr)));
        [Efield,~,~] = real_time_gpu(Q,Ax,params{1},params{2},params{3},params{4},params{5},...
            params{6},params{7},params{8},params{9},params{10},params{11},params{12},params{13},params{14},params{15},...
            params{16},params{17},params{18},params{19},params{20},params{21},params{22},params{23},params{24},...
            Tr,Tri,params{29},0);
    else
        [Efield,~] = real_time_cpu(Q,Ax,params{1},params{2},params{3},params{4},params{8},params{17},params{18},...
            params{19},params{20},params{21},params{22},Transformation,params{29});
    end
    Efield = Efield.*multiplying_factor;
end