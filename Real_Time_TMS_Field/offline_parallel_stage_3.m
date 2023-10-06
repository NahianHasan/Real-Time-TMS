function [] = offline_parallel_stage_3(NModes,FEMORD,output_folder)
    addpath(fullfile('.','matlab'))
    pathparts = strsplit(output_folder,filesep);
    subject_folder = pathparts{end};
    start_time = tic;
    for ix=1:NModes
        if exist(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(ix),'.mat']),'file')
            variableInfo = who('-file', save_file);
            status=ismember('Qi', variableInfo);
        else
            status = 0;
        end
        if status==0
            if ix>1
                disp('It seems some of the modes were deleted by someone unexpectedly.')
            end
            break;
        end
    end
    if ~status       
        F = matfile(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],[subject_folder,'_FEM_',num2str(FEMORD),'.mat']));
        volume = F.volume;
        volume = repmat(sqrt(volume),3,1);
        volume = reshape(volume,[],1);
        %Load the initially calculated modes
        Time_2 = 0;
        for i=1:NModes
            disp(['Reading Modess = ',num2str(i),'/',num2str(NModes)]);
            Q(:,i) = (load(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['coil2roi_',num2str(i),'.mat'])).Efield(:)).*volume;
            Time_2 = Time_2 + load(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['coil2roi_',num2str(i),'.mat'])).Time_2;
        end
        Time_2 = Time_2/NModes;
        %Take the QR decomposition of the mode matrix save the modes in the order of
        % highest to lowest singular values 
        [Q,R]=qr(Q,0);
        [U,~,~] = svd(R,'econ');
        U = Q*U;
        for jx=1:NModes
            U(:,jx) = U(:,jx)./volume;
        end
        [U,~] = qr(U,0);
        for ix=1:NModes
            Qi = U(:,ix);
            save(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(ix),'.mat']),'Qi','-v7.3');
        end
        %Delete the unnecessary mode files
        for i=1:NModes
            del_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['coil2roi_',num2str(i),'.mat']);
            delete(del_file)
        end
        Time_3 = toc(start_time);
        save(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],[subject_folder,'_FEM_',num2str(FEMORD),'.mat']),'Time_2','Time_3','-append')
    end
end
