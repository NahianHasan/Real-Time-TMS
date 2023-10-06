function [] = offline_parallel_stage_4(NModes,FEMORD,mode,output_folder)
    addpath(fullfile('.','matlab'))
    pathparts = strsplit(output_folder,filesep);
    subject_folder = pathparts{end};
    
    save_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Ax_',num2str(mode),'.mat']);
    if exist(save_file,'file')
        variableInfo = who('-file', save_file);
        status=ismember('Ai', variableInfo);
    else
        status = 0;
    end
    
    if ~status
        start_time = tic;
        %Load each Q mode
        F = matfile(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],[subject_folder,'_FEM_',num2str(FEMORD),'.mat']));
        disp(['ROI2COIL = ',num2str(mode),'/',num2str(NModes)]);
        Q = load(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(mode),'.mat'])).Qi;
        %Calculate surface currents
        that=reshape(Q,[3 numel(F.teid)]);
        [Js,Ks,rho,~,~]=ROItoScalp(F.te2p,F.p,F.conductivity,F.teid,that,F.rs,F.nhat,FEMORD);
        %Calculate each row basis
        Ai=cat(1,Js,Ks,rho);
        Ai = reshape(Ai,[],1);
        if ~exist(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]),'dir')
            mkdir(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]))
        end
        Time_4 = toc(start_time);
        save(save_file,'Ai','Time_4','-v7.3');
    else
        disp(['Mode (Ai) = ',num2str(mode),' is alraedy calculated'])
        return
    end
end