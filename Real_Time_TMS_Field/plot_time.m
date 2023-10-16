function [] = plot_time(current_axis,grid_spacing,Mode_Array,FEMORD,hardware,output_folders,...
                            num_simulations,subject_folders,save_dir,coil_models,plot_type,data_type)
    [Times] = collect_time(grid_spacing,Mode_Array,FEMORD,output_folders,subject_folders,save_dir,...
                                                                hardware,num_simulations,coil_models);
    plot_data(current_axis,Times,plot_type,Mode_Array,save_folder,hardware,data_type,grid_spacing)
end

function [] = plot_data(current_axis,Times,plot_type,Mode_Array,save_dir,hardware,name,grid_spacing)
    if strcmpi(name,'preprocessing_time')
        label = 'Pre-processing';
        legend_labels = {'Data Preparation','White Noise Current and Field Generation',...
            'Orthonormal Mode Function Generation',"Huygens' Surface Current Generation",'Total TIme'};
    elseif strcmpi(name,'real_time')
        label = 'Reconstruction';
        legend_labels = {"Huygens' Surface Transformation",'Primary Field Interpolation',...
            'Mode Coefficient Calculation','TMS Field Calculation','Total Time'};
    end
    if strcmpi(plot_type,'patch')
        plot_patch(current_axis,Mode_Array,Times,name,hardware,label)
    elseif strcmpi(plot_type,'mean')
        plot_mean(current_axis,Mode_Array,Times,name,hardware,label)
    elseif strcmpi(plot_type,'range')
        plot_range(current_axis,Mode_Array,Times,name,hardware,label);
    else
        disp('Please specify the plot type as either of patch, mean or, range.')
        return;
    end
    xlim([0, length(Mode_Array)+1]);
    xlabel(current_axis,'Number of Modes')
    xticks(current_axis,1:1:max(xlim))
    xticklabels(current_axis,{Mode_Array})
    xtickangle(current_axis,90);
    LGD = legend(current_axis,legend_labels);
    LGD.Location = 'northwest';
    grid(current_axis,'on');box(current_axis,'on')
    set(current_axis,'FontSize',12,'fontweight','bold','linewidth',2);
    name_2 = [label,' Time_',hardware,'_',num2str(grid_spacing)];
    exportgraphics(current_axis,fullfile(save_dir,[name_2,'_',plot_type,'.png']),...
        'ContentType','vector','Resolution',300,'BackgroundColor','w')
    saveas(current_axis,fullfile(save_dir,[name_2,'_',plot_type,'.fig']))
end

function [] = plot_patch(current_axis,Mode_Array,Times,name,hardware,label)
    X = 1:1:length(Mode_Array);
    X1 = [X, fliplr(X)];
    num_modes = length(Mode_Array);
    color_array = {'r','b','k','c','m','g'};%equal to at least the number of coil models
    marker_array = {'*','o','s','d','p','h'};%equal to at least the number of coil models
    for step=1:5
        if strcmpi(name,'real_time') 
            inBetween = [min(Times.min_real_time(:,1:num_modes,:,step),[],1), max(Times.max_real_time(:,1:num_modes,:,step),[],1)];
            fill(current_axis,X1, inBetween, color_array{cx},'HandleVisibility','off','FaceAlpha',0.2);hold on
            plot(current_axis,1:num_modes,squeeze(mean(Times.average_realtime_time(:,1:num_modes,:,step),1)),...
                ['-',marker_array{step},color_array{step}],'MarkerSize',10,'MarkerEdgeColor','k',...
                'MarkerFaceColor',color_array{step},'LineWidth',2,'DisplayName','');
            if strcmpi(hardware,'cpu')   
                ylim(current_axis,[0 3000]);
            elseif strcmpi(hardware,'gpu')
                ylim(current_axis,[0 4]);
            end
            ylabel(current_axis,[upper(hardware),' ',label,' ','Time (ms)']);
        elseif strcmpi(name,'preprocessing_time')
            inBetween = [min(Times.min_setup_time(:,1:num_modes,:,step),[],1), max(Times.max_setup_time(:,1:num_modes,:,step),[],1)];
            fill(current_axis,X1, inBetween, color_array{cx},'HandleVisibility','off','FaceAlpha',0.2);hold on
            plot(current_axis,1:length(Mode_Array),squeeze(Times.average_setup_time(1,:,:,step)),...
                ['-',marker_array{step},color_array{step}],'MarkerSize',10,'MarkerEdgeColor','k',...
                'MarkerFaceColor',color_array{step},'LineWidth',2,'DisplayName','');
            ylim(current_axis,[0 70]);
            ylabel(current_axis,[upper(hardware),' ',label,' ','Time (hr)']);
        end
    end
end


function [] = plot_range(current_axis,Mode_Array,Times,name,hardware,label)
    color_array = {'r','b','k','c','m','g'};
    marker_array = {'*','o','s','d','p','h'};
    num_modes = length(Mode_Array);
    hold(current_axis,'on');
    for step=1:5
        if strcmpi(name,'real_time') 
            errorbar(current_axis,1:num_modes,squeeze(mean(Times.average_realtime_time(:,1:num_modes,:,step),1)),...
            squeeze(mean(Times.average_realtime_time(:,1:num_modes,:,step),1))-squeeze(mean(Times.min_real_time(:,1:num_modes,:,step),1)),...
            squeeze(mean(Times.max_real_time(:,1:num_modes,:,step),1))-squeeze(mean(Times.average_realtime_time(:,1:num_modes,:,step),1)),...
            ['-',marker_array{step},color_array{step}],'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',color_array{step},...
            'LineWidth',2,'DisplayName',''); hold on
            plot(current_axis,1:num_modes,squeeze(mean(Times.average_realtime_time(:,1:num_modes,:,step),1)),...
                ['-',marker_array{step},color_array{step}],'MarkerSize',10,'MarkerEdgeColor','k',...
                'MarkerFaceColor',color_array{step},'LineWidth',2,'DisplayName','');
            if strcmpi(hardware,'cpu')   
                ylim(current_axis,[0 3000]);
            elseif strcmpi(hardware,'gpu')
                ylim(current_axis,[0 4]);
            end
            ylabel(current_axis,[upper(hardware),' ',label,' ','Time (ms)']);
        elseif strcmpi(name,'preprocessing_time')
            errorbar(current_axis,1:length(Mode_Array),squeeze(Times.average_setup_time(1,1:num_modes,:,step)),...
            squeeze(Times.average_setup_time(1,1:num_modes,:,step))-squeeze(Times.min_setup_time(1,1:num_modes,:,step)),...
            squeeze(Times.max_setup_time(1,1:num_modes,:,step))-squeeze(Times.average_setup_time(1,1:num_modes,:,step)),...
            ['-',marker_array{step},color_array{step}],'MarkerSize',10,'MarkerEdgeColor','k',...
            'MarkerFaceColor',color_array{step},'LineWidth',2,'DisplayName',''); hold on
            plot(current_axis,1:length(Mode_Array),squeeze(Times.average_setup_time(1,:,:,step)),...
                ['-',marker_array{step},color_array{step}],'MarkerSize',10,'MarkerEdgeColor','k',...
                'MarkerFaceColor',color_array{step},'LineWidth',2,'DisplayName','');
            ylim(current_axis,[0 70]);
            ylabel(current_axis,[upper(hardware),' ',label,' ','Time (hr)']);
        end
    end
end

function [] = plot_mean(current_axis,Mode_Array,Times,name,hardware,label)
    color_array = {'r','b','k','c','m','g'};
    marker_array = {'*','o','s','d','p','h'};
    num_modes = length(Mode_Array);
    hold(current_axis,'on');
    for step=1:5
        if strcmpi(name,'real_time') 
            plot(current_axis,1:num_modes,squeeze(mean(Times.average_realtime_time(:,1:num_modes,:,step),1)),...
                ['-',marker_array{step},color_array{step}],'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',...
                color_array{step},'LineWidth',2,'DisplayName',''); hold on
            if strcmpi(hardware,'cpu')   
                ylim(current_axis,[0 3000]);
            elseif strcmpi(hardware,'gpu')
                ylim(current_axis,[0 4]);
            end
            ylabel(current_axis,[upper(hardware),' ',label,' ','Time (ms)']);
        elseif strcmpi(name,'preprocessing_time')
            plot(current_axis,1:length(Mode_Array),squeeze(Times.average_setup_time(1,:,:,step)),...
                ['-',marker_array{step},color_array{step}],'MarkerSize',10,'MarkerEdgeColor','k',...
                'MarkerFaceColor',color_array{step},'LineWidth',2,'DisplayName',''); hold on
            ylim(current_axis,[0 70]);
            ylabel(current_axis,[upper(hardware),' ',label,' ','Time (hr)']);
        end
    end
end
