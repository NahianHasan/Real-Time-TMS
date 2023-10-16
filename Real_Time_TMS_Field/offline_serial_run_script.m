function [] = offline_serial_run_script(msh_file_read_fcn,msh_file_read_fcn_location,NModes,msh_file,FEMORD,output_folder)
    start_time = tic;
    addpath(msh_file_read_fcn_location);
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
    [rs,nhat,tri_hg,extruded_scalp_points,area_tri]=generateextrudedmesh(p,te2p);
    %Generate random mode weight matrix for the basis fields
    Wn=randn([numel(rs) NModes]);
    %calculate the volume of tetrahedrons
    nt=size(te2p,2);
    edges=reshape(p(:,reshape(te2p(2:4,:),nt*3,1)),3,3,nt)-repmat(reshape(p(:,te2p(1,:)),3,1,nt),[1,3,1]);
    volume=squeeze(1/6*sqrt(sum(dot(cross(edges(:,1,:),edges(:,2,:),1),edges(:,3,:),1).^2,1)));
    volume = volume(teid);
    volume = repmat(sqrt(volume),3,1);
    volume = reshape(volume,[],1);
    
    Time_1 = toc(start_time);
    save(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],[subject_folder,'_FEM_',num2str(FEMORD),'.mat']),'Time_1','te2p','scalp_tri','scalp_points','teid','p','reg','conductivity','NModes','FEMORD','ro','rs','nhat','area_tri','tri_hg','extruded_scalp_points','Wn','volume','-v7.3');

    %%%%%%%%%%%%%%%%%%  Stage 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ix=1:NModes
        save_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(ix),'.mat']);
        if exist(save_file,'file')
            variableInfo = who('-file', save_file);
            status=ismember('Qi', variableInfo);
        else
            status = 0;
        end
        if ~status
            if ix>1
                disp('It seems some of the modes were deleted by someone unexpectedly. To correctly calculate the modes, starting again.')
            end
            for jx=1:NModes
                del_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Ax_',num2str(jx),'.mat']);
                if exist(del_file,'file')
                    delete(del_file)
                end
            end
            break;
        end
    end
    if ~status
        Time_2e = [];
        for ix=1:NModes
            save_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['coil2roi_',num2str(ix),'.mat']);
            start_time = tic();
            disp(['Coil2ROI = ',num2str(ix),'/',num2str(NModes)]);
            if exist(save_file,'file')
                variableInfo = who('-file', save_file);
                status=ismember('Efield', variableInfo);
            else
                status = 0;
            end
            if ~status
                %define the huygens surface dipole with random weight
                ks=reshape(Wn(:,ix),size(rs));
                %calculate the E-field mode in the ROI
                Efield=runcodeks(te2p,p,conductivity,rs,ks,teid,FEMORD);
                Time_2 = toc(start_time);
                save(save_file,'Efield','Time_2');
                Q(:,ix) = Efield(:).*volume;
            else
                Q(:,ix) = (load(save_file).Efield(:)).*volume;
                Time_2 = load(save_file).Time_2;
            end
            Time_2e = [Time_2e, Time_2];
        end
        
   
        %%%%%%%%%%%%%%%%%%  Stage 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~exist(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]),'dir')
            mkdir(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]))
        end
        start_time = tic();
        %Take the QR decomposition of the mode matrix save the modes in the order of
        % highest to lowest singular values 
        [Q,R]=qr(Q,0);
        [U,~,~] = svd(R,'econ');
        U = Q*U;
        for jx=1:NModes
            U(:,jx) = U(:,jx)./volume;
        end
        [U,~] = qr(U,0);
        Q=[];
        for ix=1:NModes
            Qi = U(:,ix);
            Q(:,ix) = Qi;
            Time_2 = Time_2e(ix);
            save(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(ix),'.mat']),'Qi','Time_2','-v7.3');
        end
        Time_3 = toc(start_time);
        save(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],[subject_folder,'_FEM_',num2str(FEMORD),'.mat']),'Time_3','-append')
        %Delete the unnecessary mode files
        for ix=1:NModes
            del_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['coil2roi_',num2str(ix),'.mat']);
            delete(del_file)
        end
    else
        for ix=1:NModes
            Q(:,ix) = load(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(ix),'.mat'])).Qi;
        end
    end
    
    %%%%%%%%%%%%%%%%%%  Stage 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ix=1:NModes
        save_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Ax_',num2str(ix),'.mat']);
        if exist(save_file,'file')
            variableInfo = who('-file', save_file);
            status=ismember('Ai', variableInfo);
        else
            status = 0;
        end
        if ~status
            disp(['ROI2COIL = ',num2str(ix),'/',num2str(NModes)]);
            start_time = tic();
            %Calculate surface currents
            that=reshape(Q(:,ix),[3 numel(teid)]);
            [Js,Ks,rho,~,~]=ROItoScalp(te2p,p,conductivity,teid,that,rs,nhat,FEMORD);
            %Calculate each row basis
            Ai=cat(1,Js,Ks,rho);
            Ai = reshape(Ai,[],1);
            Time_4 = toc(start_time);
            save(save_file,'Ai','Time_4','-v7.3');
        else
            continue;
        end
    end
    disp([newline,'Mode generration complete.',newline])
end
