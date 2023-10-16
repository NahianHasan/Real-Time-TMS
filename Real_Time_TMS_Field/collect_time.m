function [Times] = collect_time(grid_spacing,Mode_Array,FEMORD,output_folders,subject_folders,save_folder,...
                                  hardware,num_simulations,coil_models)
    percentile= 100;
    Offline_Setup_Times = zeros([1,length(Mode_Array),length(output_folders),1,4]);
    Realtime_Times = zeros([length(coil_models),length(Mode_Array),length(output_folders),num_simulations,4]);
    Realtime_Setup_Times = zeros([length(coil_models),length(Mode_Array),length(output_folders),num_simulations,1]);
    Realtime_Communication_Times = zeros([length(coil_models),length(Mode_Array),length(output_folders),num_simulations,1]);
    %collect real-time times
    for cx=1:length(coil_models)
        coil_model = coil_models{cx};
        for jx=1:length(output_folders)
            output_folder = output_folders{jx};
            for ix=1:length(Mode_Array)
                NModes = Mode_Array(ix);
                W = load(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Timing_grid_',num2str(grid_spacing),'_Modes_',num2str(NModes),'_coil_',coil_model,'_',lower(hardware),'.mat']));
                temp = squeeze(W.Time);
                Realtime_Times(cx,ix,jx,:,:) = [temp(2:end,1),temp(2:end,2),temp(2:end,3),temp(2:end,4)];
                Realtime_Setup_Times(cx,ix,jx,:,:) = W.data_set_up_time;
                Realtime_Communication_Times(cx,ix,jx,:,:) = W.communication_time;
            end
        end
    end
    if ~exist(fullfile(sae_dir,'Collected_Setup_Time_data.mat'),'file')
        for jx=1:length(output_folders)
            output_folder = output_folders{jx};
            subject_folder = subject_folders{jx};
            for ix=1:length(Mode_Array)
                NModes = Mode_Array(ix);
                W = matfile(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],[subject_folder,'_FEM_',num2str(FEMORD),'.mat']));
                T4 = 0;
                T2 = 0;
                for kx=1:NModes
                    M = matfile(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(kx),'.mat']));
                    T2 = T2 + M.Time_2;
                    M = matfile(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Ax_',num2str(kx),'.mat']));
                    T4 = T4 + M.Time_4;
                end
                Tm = [W.Time_1,T2,W.Time_3,T4];
                Offline_Setup_Times(1,ix,jx,:,:) = Tm;
                disp(['Subject = ',subject_folder,', Mode = ',num2str(NModes)])
            end
            save(fullfile(save_folder,'Collected_Setup_Time_data.mat'),'Offline_Setup_Times','-v7.3')
        end
    else
        Setup_Times = load(fullfile(save_folder,'Collected_Setup_Time_data.mat')).Offline_Setup_Times;
    end

    Realtime_Times = Realtime_Times.*1000;%convert time to milliseconds
    Offline_Setup_Times = Setup_Times./3600;%convert time to hours

    Realtime_Times = reshape(Realtime_Times,size(Realtime_Times,1),size(Realtime_Times,2),[],size(Realtime_Times,5));
    Total_Real_Time = sum(Realtime_Times,4);
    Realtime_Times = cat(4,Realtime_Times,Total_Real_Time);
    average_realtime_time = mean(Realtime_Times,3);
    max_real_time = prctile(Realtime_Times,percentile,3);
    min_real_time = min(Realtime_Times,[],3);
    Offline_Setup_Times = reshape(Offline_Setup_Times,size(Offline_Setup_Times,1),size(Offline_Setup_Times,2),[],size(Offline_Setup_Times,5));
    Total_setup_Time = sum(Offline_Setup_Times,4);
    Offline_Setup_Times = cat(4,Offline_Setup_Times,Total_setup_Time);
    average_setup_time = mean(Offline_Setup_Times,3);
    max_setup_time = max(Offline_Setup_Times,[],3);
    min_setup_time = min(Offline_Setup_Times,[],3);
    Times = struct();
    Times.Offline_Setup_Times = Offline_Setup_Times;
    Times.Realtime_Times = Realtime_Times;
    Times.average_realtime_time = average_realtime_time;
    Times.average_setup_time = average_setup_time;
    Times.max_real_time = max_real_time;
    Times.min_real_time = min_real_time;
    Times.max_setup_time = max_setup_time;
    Times.min_setup_time = min_setup_time;
    Times.Realtime_Setup_Times = Realtime_Setup_Times;
    Times.Realtime_Communication_Times = Realtime_Communication_Times;
end

