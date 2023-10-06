function [rcoil,kcoil,coil_model,th_hair,coil_tri,base_positions] = load_coil_model(real_time_code_path,coil_model_file,save_model)
    clear rcoil kcoil
    M = load(coil_model_file);
    t = strsplit(coil_model_file,filesep);
    coil_model = t{end}(1:end-4);
    if strcmpi(coil_model,'coil_fig-8')
        rcoil = M.rcoil;%shape = N*3
        kcoil = M.jcoil;%shape = N*3
        coil_tri = boundary(rcoil(:,1),rcoil(:,2),rcoil(:,3));
        base_positions = rcoil;
        th_hair = 0.005;
    elseif strcmpi(coil_model,'coil_DB80-2')
        base_positions = M.P;
        rcoil = M.rv;%shape = N*3
        kcoil = M.Jdc;%shape = N*3 
        coil_tri = M.t;
        th_hair=0.02;
    elseif strcmpi(coil_model,'coil_Cool40')
        rcoil = M.rv;%shape = N*3
        kcoil = M.Jdc;%shape = N*3 
        coil_tri = M.t;
        base_positions = rcoil;
        th_hair=0.015;
    elseif strcmpi(coil_model,'coil_CB60')
        rcoil = M.rv;%shape = N*3
        kcoil = M.Jdc;%shape = N*3 
        coil_tri = M.t;
        base_positions = rcoil;
        th_hair=0.005;
    end
    if save_model
        if ~exist(fullfile(real_time_code_path,'coil_data.mat'),'file')
            save(fullfile(real_time_code_path,'coil_data.mat'),'coil_model','rcoil','kcoil','coil_tri','th_hair','-v7.3')
        end
    end
end