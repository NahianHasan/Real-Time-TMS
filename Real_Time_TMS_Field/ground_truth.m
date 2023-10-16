function [E_org,Time_gt,Transformation] = ground_truth(real_time_code_path,mode,msh_file,msh_file_read_fcn,...
            msh_file_read_fcn_location,m2m_dir,coil_model_file,FEMORD,mapping_surface,output_folder,varargin)    
    multiplying_factor = 6.6E7;
    addpath(msh_file_read_fcn_location);
    addpath(fullfile(real_time_code_path,'matlab'));
    if exist(fullfile(m2m_dir,'mri2mesh_log.html'),'file')
        model_creation_tool = 'mri2mesh';
    elseif exist(fullfile(m2m_dir,'headreco_log.html'),'file')
        model_creation_tool = 'headreco';
    end
    
    [rcoil,kcoil,coil_model,~,~,~] = load_coil_model(real_time_code_path,coil_model_file,0);
    [p,te2p,conductivity,~] = load_msh_data(msh_file,msh_file_read_fcn);
    
    if size(varargin)==1
        Transformations = varargin{1};
    else
        data_path = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Random_Coil_Placement_IDs_',coil_model,'.mat']);
        Transformations = load(data_path).Transformations;
    end
    
    %select only those tetrahedrons within the GM/WM - Offline ROI
    teid=1:numel(te2p)/4;
    teid=teid(conductivity==0.2750 | conductivity==0.1260);
    
    %find the teids corresponding to GM surface or GM/WM middle layer
    [GM_teid,~] = get_mapping_surface(msh_file,msh_file_read_fcn,mapping_surface,m2m_dir,model_creation_tool);
    
    if ~exist(fullfile(output_folder,['FEM_',num2str(FEMORD)],['GT_E_Fields_',coil_model]),'dir')
        mkdir(fullfile(output_folder,['FEM_',num2str(FEMORD)],['GT_E_Fields_',coil_model]))
    end
    
    if size(Transformations,3)==1
        ID = 1;
    else
        ID = mode;
    end
    Transformation = Transformations(:,:,ID);
    start_time = tic();
    rs = rcoil';
    ks = kcoil';
    rs(4,:)=1;
    rs2 = Transformation(1:3,1:4)*rs;
    ks2=Transformation(1:3,1:3)*ks;
    [E_org,~,~]=runcode(te2p,p,conductivity,rs2,ks2,teid,FEMORD);
    E_org = E_org.*multiplying_factor;
    if size(E_org,1)<3
        E_org = reshape(E_org,3,[]);
    end
    E_org = E_org(:,GM_teid);
    Time_gt = toc(start_time);
    save(fullfile(output_folder,['FEM_',num2str(FEMORD)],['GT_E_Fields_',coil_model],['E_org_',num2str(mode),'.mat']),'E_org','Time_gt','Transformation','-v7.3')
end
