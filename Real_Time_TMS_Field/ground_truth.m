function [E_org,Time_gt,Transformation] = ground_truth(real_time_code_path,mode,msh_file,msh_file_read_fcn,msh_file_read_fcn_location,m2m_dir,coil_model_file,FEMORD,mapping_surface,output_folder,varargin)    
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
    if strcmp(mapping_surface, 'GM')% GM/WM surface
        xx = [];
        xx(conductivity==.1260)=1; % White
        xx(conductivity==.2750)=1; % Grey
        [~,GM_teid]=surftri(p',te2p(:,xx(:)==1)'); % GM/WM surface
        %[GM_tri,GM_teid]=surftri(p',te2p(:,conductivity==0.2750)'); 
    elseif strcmp(mapping_surface, 'GMM')% GM middle surface
        [~,GM_msh] = load_GM_mid_Layer(m2m_dir,model_creation_tool);
        v1 = GM_msh.nodes(GM_msh.triangles(:,2),:)-GM_msh.nodes(GM_msh.triangles(:,1),:);
        v2 = GM_msh.nodes(GM_msh.triangles(:,3),:)-GM_msh.nodes(GM_msh.triangles(:,1),:);
        SA = sum(0.5*vecnorm(cross(v1,v2),2,2),'all')/1E6;
        warning('off');TR = triangulation(te2p(:,teid)',p');warning('on')
        GM_teid = pointLocation(TR,GM_msh.nodes/1000);
        GM_teid(isnan(GM_teid)) = 1;
        clear TR
    end
    
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
