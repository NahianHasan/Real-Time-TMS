function [] = offline_parallel_stage_2(NModes,FEMORD,mode,output_folder)
    addpath(fullfile('.','matlab'))
    pathparts = strsplit(output_folder,filesep);
    subject_folder = pathparts{end};
    
    save_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['coil2roi_',num2str(mode),'.mat']);
    if exist(save_file,'file')
        variableInfo = who('-file', save_file);
        status=ismember('Efield', variableInfo);
    else
        if exist(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(mode),'.mat']),'file')
            variableInfo = who('-file', fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(mode),'.mat']));
            status=ismember('Qi', variableInfo);
        else
            status = 0;
        end
    end
    if ~status
        start_time = tic();
        %Load subject data
        F = matfile(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],[subject_folder,'_FEM_',num2str(FEMORD),'.mat']));
        disp(['Coil2ROI = ',num2str(mode),'/',num2str(NModes)]);
        rs = F.rs; Wn = F.Wn; p = F.p; te2p = F.te2p; teid = F.teid; conductivity = F.conductivity;
        %define the huygens surface dipole with random weight
        ks=reshape(Wn(:,mode),size(rs));
        %calculate the E-field mode in the ROI
        Efield=runcodeks(te2p,p,conductivity,rs,ks,teid,FEMORD);
        Time_2 = toc(start_time);
        if ~exist(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]),'dir')
            mkdir(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]))
        end
        save(save_file,'Efield','Time_2');
    else
        disp(['Mode (Q) = ',num2str(mode),' is alraedy calculated'])
        return
    end
end