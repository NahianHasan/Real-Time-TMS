function [cluster_params] = collect_cluster_parameters(file_name)
    cluster_params = {};
    cnt = 1;
    F=fopen(file_name,'r');
    tline = fgetl(F);
    while ischar(tline)
        L=tline;
        t=strsplit(L,',');
        cluster_params{cnt} = t{2};
        cnt = cnt + 1;
        tline = fgetl(F);
    end
end

