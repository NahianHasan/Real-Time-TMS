function [Q,Ax,params] = real_time_stage_setup(msh_file,msh_file_read_fcn,real_time_code_path,...
        NModes,m2m_dir,FEMORD,grid_spacing,mapping_surface,output_folder,...
                        subject_folder,coil_model_file)
    start_t = tic();
    if exist(fullfile(m2m_dir,'mri2mesh_log.html'),'file')
        model_creation_tool = 'mri2mesh';
    elseif exist(fullfile(m2m_dir,'headreco_log.html'),'file')
        model_creation_tool = 'headreco';
    end
    
    f = waitbar(0,'Loading Modes...');
    
    [~,~,coil_model,~,~,~] = load_coil_model(real_time_code_path,coil_model_file,0);
    grid_fields_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['grid_fields_',num2str(grid_spacing),'_',coil_model,'.mat']);
    if ~exist(grid_fields_file,'file')
        disp(['Primary fields are not calculated, generating primary fields for the coil model ', coil_model])
        primary_field_generation(msh_file,msh_file_read_fcn,real_time_code_path,coil_model_file,grid_spacing,FEMORD,output_folder);
    end
    %find whether there is gpu in the system or not
    try
        g = gpuDevice;
        if g.DeviceAvailable
            gpu_solution = 1;
            reset(g);
        else
            g = [];
            gpu_solution = 0;
        end
    catch
        g = [];
        gpu_solution = 0;
    end
    interp_method = 'linear';

    %%%%%%%%%% Find the mapping surface tetra IDs %%%%%%%%%%%%%%%%%%%%
    G1 = matfile(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],[subject_folder,'_FEM_',num2str(FEMORD),'.mat']));
    rs = G1.rs;
    [GM_teid,GM_msh] = get_mapping_surface(msh_file,msh_file_read_fcn,mapping_surface,m2m_dir,model_creation_tool);
    %%%%%%%%%% Load offline modes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q = zeros(3*numel(GM_teid),NModes);
    disp('Loading Mode Functions');
    for ix=1:NModes
        Qt = load(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(ix),'.mat'])).Qi;
        Qt = reshape(Qt,3,[]);
        Qt = Qt(:,GM_teid);
        Q(:,ix) = reshape(Qt,[],1);
        waitbar(ix/NModes,f,sprintf('Number of E-Field Modes = %4.0f/%4.0f',ix,NModes))
    end
    disp(['Loaded modes = ',num2str(NModes)]);
    disp('Loading Surface Current Modes');
    Ax_temp = zeros(7*numel(rs)/3,NModes);
    for ix=1:NModes
        At = load(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Ax_',num2str(ix),'.mat'])).Ai;
        Ax_temp(:,ix) = At';
        waitbar(ix/NModes,f,sprintf('Number of Current Modes = %4.0f/%4.0f',ix,NModes))
    end
    disp(['Loaded Surface Current Modes = ',num2str(NModes)]);
    Ax_temp = reshape(Ax_temp(:,1:NModes),[7,numel(rs)/3,NModes]);
    close(f)

    %%%%%%%%%% Load pre-computed grid fields %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G = load(grid_fields_file);
    X=G.X;Y=G.Y;Z=G.Z;
    Eout = G.Eout;
    Hout = G.Hout;
    size_X = size(X);


    %%%%%%%%%% Load data to GPU if possible  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    huygens_surf_points = rs;
    huygens_surf_points(4,:)=1;
    huygens_surf_points = single(huygens_surf_points);
    huygens_surf_points_rot_x = single(zeros(1,size(huygens_surf_points,2)));
    huygens_surf_points_rot_y = single(zeros(1,size(huygens_surf_points,2)));
    huygens_surf_points_rot_z = single(zeros(1,size(huygens_surf_points,2)));
    NModes = single(NModes);
    Exv = single(reshape(Eout(1,:),size_X));
    Eyv = single(reshape(Eout(2,:),size_X));
    Ezv = single(reshape(Eout(3,:),size_X));
    Hxv = single(reshape(Hout(1,:),size_X));
    Hyv = single(reshape(Hout(2,:),size_X));
    Hzv = single(reshape(Hout(3,:),size_X));
    X = single(X);
    Y = single(Y);
    Z = single(Z);
    Tr = single(zeros(4,4));
    Tri = single(zeros(4,4));
    SEx = single(zeros(size(huygens_surf_points_rot_x,2),1));
    SEy = single(zeros(size(huygens_surf_points_rot_x,2),1));
    SEz = single(zeros(size(huygens_surf_points_rot_x,2),1));
    SHx = single(zeros(size(huygens_surf_points_rot_x,2),1));
    SHy = single(zeros(size(huygens_surf_points_rot_x,2),1));
    SHz = single(zeros(size(huygens_surf_points_rot_x,2),1));
    Q = single(Q(:,1:NModes));
    hardware = 'cpu';
    if gpu_solution
        hardware = 'gpu';
        %dynamically define how many sub arrays you need for efficient GPU
        %computation. We consider 120,000*100 to be the maximum size of a
        %sub_matrix Axi
        Nh = size(huygens_surf_points_rot_x,2);%number of huygens surf points
        k = 12000000;
        gx = floor(k/Nh);
        Ax = {};
        Fields = {};
        coeff = {};
        count = 1;
        ix = 1;
        while(1)
            if (count+gx-1) <= NModes
                Ax{ix} = single(Ax_temp(1:6,:,count:count+gx-1));
                Fields{ix} = single(zeros(6,size(huygens_surf_points_rot_x,2),gx));
                coeff{ix} = single(zeros(gx,1));
            else
                Ax{ix} = single(Ax_temp(1:6,:,count:NModes));
                Fields{ix} = single(zeros(6,size(huygens_surf_points_rot_x,2),NModes-count+1));
                coeff{ix} = single(zeros(NModes-count+1,1));
            end
            count = count + gx;
            ix = ix+1;
            if count > NModes
                break;
            end
        end
    else
        Ax = Ax_temp;
        Fields = {};
        coeff = {};
    end
    clear Ax_temp;
    Fields_temp = single(zeros(6,size(huygens_surf_points_rot_x,2)));
    Efield = single(zeros(size(Q,1),1));
    data_set_up_time = toc(start_t);
    if gpu_solution
        start_t = tic();
        X = gpuArray(X);
        Y = gpuArray(Y);
        Z = gpuArray(Z);
        Exv = gpuArray(Exv);
        Eyv = gpuArray(Eyv);
        Ezv = gpuArray(Ezv);
        Hxv = gpuArray(Hxv);
        Hyv = gpuArray(Hyv);
        Hzv = gpuArray(Hzv);
        NModes = gpuArray(NModes);
        for ix=1:length(coeff)
            coeff{ix} = gpuArray(coeff{ix});
        end
        Tr = gpuArray(Tr);
        Tri = gpuArray(Tri);
        
        huygens_surf_points = gpuArray(huygens_surf_points);
        huygens_surf_points_rot_x = gpuArray(huygens_surf_points_rot_x);
        huygens_surf_points_rot_y = gpuArray(huygens_surf_points_rot_y);
        huygens_surf_points_rot_z = gpuArray(huygens_surf_points_rot_z);
        SEx = gpuArray(SEx);
        SEy = gpuArray(SEy);
        SEz = gpuArray(SEz);
        SHz = gpuArray(SHx);
        SHy = gpuArray(SHy);
        SHx = gpuArray(SHz);
        
        Q = gpuArray(Q);
        for ix=1:length(coeff)
            Ax{ix} = gpuArray(Ax{ix});
        end
        Fields_temp = gpuArray(Fields_temp);
        for ix=1:length(coeff)
            Fields{ix} = gpuArray(Fields{ix});
        end
        Efield = gpuArray(Efield);
        communication_time = toc(start_t);
    else
        communication_time = 0;
    end
    data_set_up_time = data_set_up_time + communication_time;
    params = {};
    params{1} = X;
    params{2} = Y;
    params{3} = Z;
    params{4} = huygens_surf_points;
    params{5} = huygens_surf_points_rot_x;
    params{6} = huygens_surf_points_rot_y;
    params{7} = huygens_surf_points_rot_z;
    params{8} = interp_method;
    params{9} = Fields;
    params{10} = Fields_temp;
    params{11} = SEx;
    params{12} = SEy;
    params{13} = SEz;
    params{14} = SHx;
    params{15} = SHy;
    params{16} = SHz;
    params{17} = Exv;
    params{18} = Eyv;
    params{19} = Ezv;
    params{20} = Hxv;
    params{21} = Hyv;
    params{22} = Hzv;
    params{23} = coeff;
    params{24} = Efield;
    params{25} = Tr;
    params{26} = Tri;
    params{27} = GM_teid;
    params{28} = GM_msh;
    params{29} = NModes;
    params{30} = gpu_solution;
    params{31} = g;
    params{32} = data_set_up_time;
    params{33} = communication_time;
    params{34} = coil_model;
    params{35} = hardware;
    %a dummy call to the gpu functions for initializing it
    [~] = real_time_field_calculation_single_placements(Q,Ax,params,rand([4,4]));
end
