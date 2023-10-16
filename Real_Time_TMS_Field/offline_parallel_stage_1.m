function [] = offline_parallel_stage_1(NModes,msh_file,msh_file_read_fcn,FEMORD,output_folder,msh_file_read_fcn_location)
    start_time = tic;
    addpath(msh_file_read_fcn_location);
    addpath(fullfile('.','matlab'));
    pathparts = strsplit(output_folder,filesep);
    subject_folder = pathparts{end};

    if ~exist(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]),'dir')
        mkdir(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]))
    end
    
    [p,te2p,conductivity,reg] = load_msh_data(msh_file,msh_file_read_fcn);
    
    %select only those tetrahedrons within the GM/WM - Offline ROI
    teid=1:numel(te2p)/4;
    teid=teid(conductivity==0.2750 | conductivity==0.1260);

    scalp_tri=surftri(p',te2p');%select scalp triangular facets
    scalp_points = p(:,unique(scalp_tri));%select scalp nodes
    %select the observation points as centers of ROI tetrahedra
    ro = (p(:,te2p(1,teid))+p(:,te2p(2,teid))+p(:,te2p(3,teid))+p(:,te2p(4,teid)))'./4;
    %Generate the Huygens Surface Points
    [rs,nhat,extruded_scalp_points,surf_point_ind,area_tri]=generateextrudedmesh(p,te2p);
    %Generate random mode weight matrix for the basis fields
    Wn=randn([numel(rs) NModes]);
    %calculate the volume of tetrahedrons
    nt=size(te2p,2);
    edges=reshape(p(:,reshape(te2p(2:4,:),nt*3,1)),3,3,nt)-repmat(reshape(p(:,te2p(1,:)),3,1,nt),[1,3,1]);
    volume=squeeze(1/6*sqrt(sum(dot(cross(edges(:,1,:),edges(:,2,:),1),edges(:,3,:),1).^2,1)));
    volume = volume(teid);

    Time_1 = toc(start_time);
    save(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],[subject_folder,'_FEM_',num2str(FEMORD),'.mat']),'Time_1','te2p','scalp_tri','scalp_points','teid','p','reg','conductivity','NModes','FEMORD','ro','rs','nhat','area_tri','extruded_scalp_points','surf_point_ind','Wn','volume','-v7.3');
end

