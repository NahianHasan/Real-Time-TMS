function [] = plot_memory(current_axis,grid_spacing,Mode_Array,FEMORD,output_folders,save_folder,...
                                    plot_type,coil_models)
    [Memory] = collect_memory(grid_spacing,Mode_Array,FEMORD,output_folders,coil_models);
    plot_data(current_axis,Memory,plot_type,Mode_Array,save_folder)
end

function [] = plot_data(current_axis,Memory,plot_type,Mode_Array,save_folder)
    if strcmp(plot_type,'patch')
        plot_patch(current_axis,Mode_Array,Memory)
    elseif strcmp(plot_type,'mean')
        plot_mean(current_axis,Mode_Array,Memory)
    elseif strcmp(plot_type,'range')
        plot_range(current_axis,Mode_Array,Memory)
    else
        disp('Please specify the plot type as either of patch, mean or, range.')
        return;
    end
    xlim([0, length(Mode_Array)+1]);
    xlabel(current_axis,'Number of Modes')
    xticks(current_axis,1:1:max(xlim))
    xticklabels(current_axis,{Mode_Array})
    xtickangle(current_axis,90);
    ylabel(current_axis,['Memory (GB)']);
    ylim(current_axis,[0,8]);
    LGD = legend(current_axis);
    LGD.Location = 'northwest';
    grid(current_axis,'on');box(current_axis,'on')
    set(current_axis,'FontSize',12,'fontweight','bold','linewidth',2);
    exportgraphics(current_axis,fullfile(save_folder,['Memory_',plot_type,'.png']),...
        'ContentType','vector','Resolution',300,'BackgroundColor','w')
    saveas(current_axis,fullfile(save_folder,['Memory_',plot_type,'.fig']))
end

function [] = plot_range(current_axis,Mode_Array,Memory)
    num_modes = length(Mode_Array);
    color_array = {'r','b','k','c','m','g'};
    marker_array = {'*','o','s','d','p','h'};
    hold(current_axis,'on');
    errorbar(current_axis,1:num_modes,Memory.cpu1_memory_during_real_time_in_cpu.average_memory,...
        Memory.cpu1_memory_during_real_time_in_cpu.average_memory-Memory.cpu1_memory_during_real_time_in_cpu.min_memory,...
        Memory.cpu1_memory_during_real_time_in_cpu.max_memory-Memory.cpu1_memory_during_real_time_in_cpu.average_memory,...
        ['-',marker_array{1},color_array{1}],'MarkerSize',10,'MarkerEdgeColor',color_array{1},...
        'MarkerFaceColor',color_array{1},'LineWidth',2,'DisplayName','CPU memory during real-time in CPU'); hold on
    errorbar(current_axis,1:num_modes,Memory.cpu1_memory_during_real_time_in_gpu.average_memory,...
        Memory.cpu1_memory_during_real_time_in_gpu.average_memory-Memory.cpu1_memory_during_real_time_in_gpu.min_memory,...
        Memory.cpu1_memory_during_real_time_in_gpu.max_memory-Memory.cpu1_memory_during_real_time_in_gpu.average_memory,...
        ['-',marker_array{2},color_array{2}],'MarkerSize',10,'MarkerEdgeColor',color_array{2},...
        'MarkerFaceColor',color_array{2},'LineWidth',2,'DisplayName','CPU memory during real-time in GPU'); hold on
    errorbar(current_axis,1:num_modes,Memory.gpu1_memory_during_real_time.average_memory,...
        Memory.gpu1_memory_during_real_time.average_memory-Memory.gpu1_memory_during_real_time.min_memory,...
        Memory.gpu1_memory_during_real_time.max_memory-Memory.gpu1_memory_during_real_time.average_memory,...
        ['-',marker_array{3},color_array{3}],'MarkerSize',10,'MarkerEdgeColor',color_array{3},...
        'MarkerFaceColor',color_array{3},'LineWidth',2,'DisplayName','GPU memory during real-time in GPU');
end


function [] = plot_mean(current_axis,Mode_Array,Memory)
    num_modes = length(Mode_Array);
    color_array = {'r','b','k','c','m','g'};
    marker_array = {'*','o','s','d','p','h'};
    hold(current_axis,'on');
    plot(current_axis,1:num_modes,Memory.cpu1_memory_during_real_time_in_cpu.average_memory,...
        ['-',marker_array{1},color_array{1}],'MarkerSize',10,'MarkerEdgeColor','k',...
        'MarkerFaceColor',color_array{1},'LineWidth',2,'DisplayName','CPU memory during real-time in CPU');
    plot(current_axis,1:num_modes,Memory.cpu1_memory_during_real_time_in_gpu.average_memory,...
        ['-',marker_array{2},color_array{2}],'MarkerSize',10,'MarkerEdgeColor','k',...
        'MarkerFaceColor',color_array{2},'LineWidth',2,'DisplayName','CPU memory during real-time in GPU');
    plot(current_axis,1:num_modes,Memory.gpu1_memory_during_real_time.average_memory,...
        ['-',marker_array{3},color_array{3}],'MarkerSize',10,'MarkerEdgeColor','k',...
        'MarkerFaceColor',color_array{3},'LineWidth',2,'DisplayName','GPU memory during real-time in GPU');
end

 
function [] = plot_patch(current_axis,Mode_Array,Memory)
    num_modes = length(Mode_Array);
    X = 1:1:length(Mode_Array);
    X1 = [X, fliplr(X)];
    color_array = {'r','b','k','c','m','g'};
    marker_array = {'*','o','s','d','p','h'};
    inBetween = [min(Memory.cpu1_memory_during_real_time_in_gpu.min_memory,[],1), ...
        max(fliplr(Memory.cpu1_memory_during_real_time_in_gpu.max_memory),[],1)];
    fill(current_axis,X1, inBetween, color_array{cx},'HandleVisibility','off','FaceAlpha',0.2);
    plot(current_axis,1:num_modes,Memory.cpu1_memory_during_real_time_in_gpu.average_memory,...
        ['-',marker_array{2},color_array{2}],'MarkerSize',10,'MarkerEdgeColor','k',...
        'MarkerFaceColor',color_array{2},'LineWidth',2,'DisplayName','CPU memory during real-time in GPU');

    inBetween = [min(Memory.cpu1_memory_during_real_time_in_cpu.min_memory,[],1), ...
        max(fliplr(Memory.cpu1_memory_during_real_time_in_cpu.max_memory),[],1)];
    fill(current_axis,X1, inBetween, color_array{cx},'HandleVisibility','off','FaceAlpha',0.2);
    plot(current_axis,1:num_modes,Memory.cpu1_memory_during_real_time_in_cpu.average_memory,...
        ['-',marker_array{2},color_array{2}],'MarkerSize',10,'MarkerEdgeColor','k',...
        'MarkerFaceColor',color_array{2},'LineWidth',2,'DisplayName','CPU memory during real-time in CPU');
    
    inBetween = [min(Memory.gpu1_memory_during_real_time.min_memory,[],1), ...
        max(fliplr(Memory.gpu1_memory_during_real_time.max_memory),[],1)];
    fill(current_axis,X1, inBetween, color_array{cx},'HandleVisibility','off','FaceAlpha',0.2);
    plot(current_axis,1:num_modes,Memory.gpu1_memory_during_real_time.average_memory,...
        ['-',marker_array{2},color_array{2}],'MarkerSize',10,'MarkerEdgeColor','k',...
        'MarkerFaceColor',color_array{2},'LineWidth',2,'DisplayName','GPU memory during real-time in GPU');
end
