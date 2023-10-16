function [Errors,current_axis] = plot_error(current_axis,real_time_code_path,grid_spacing,Mode_Array,FEMORD,hardware,...
                                            output_folders,save_dir,coil_model_files,plot_type,data_type)
    if ~exist(save_dir,'dir')
        mkdir(save_dir);
    end
    coil_models = {};
    for ix=1:length(coil_model_files)
        [~,~,coil_model,~,~,~] = load_coil_model(real_time_code_path,coil_model_files{ix},0);
        coil_models{ix} = coil_model;
    end
    [Errors] = collect_error(grid_spacing,Mode_Array,FEMORD,output_folders,hardware,coil_models);
    
    if strcmpi(plot_type,'patch')
        plot_patch(current_axis,Mode_Array,Errors,data_type,coil_models);
    elseif strcmpi(plot_type,'box')
        plot_box(current_axis,Mode_Array,Errors,data_type,coil_models);
    elseif strcmpi(plot_type,'range')
        plot_range(current_axis,Mode_Array,Errors,data_type,coil_models);
    elseif strcmpi(plot_type,'mean')
        plot_mean(current_axis,Mode_Array,Errors,data_type,coil_models);
    else
        disp('Please specify the plot type as either of patch, mean, box or, range.')
        return;
    end
    xlim([0, length(Mode_Array)+1]);
    ylim(current_axis,[1,10]);
    xlabel(current_axis,'Number of Modes','fontweight','bold')
    xticks(current_axis,1:1:max(xlim))
    xticklabels(current_axis,{Mode_Array})
    xtickangle(current_axis,90);
    yticks(current_axis,[1,2,3,5,10,100]);
    yticklabels(current_axis,{[1,2,3,5,10,100]})
    grid(current_axis,'on')
    LGD = legend(current_axis);
    set(current_axis,'FontSize',12,'fontweight','bold','linewidth',2)
    exportgraphics(current_axis,fullfile(save_dir,['Errors_',num2str(grid_spacing),'_',data_type,'_',plot_type,'.png']),...
        'ContentType','vector','Resolution',300,'BackgroundColor','w')
end

function [] = plot_patch(current_axis,Mode_Array,Errors,coil_models)
    X = 1:1:length(Mode_Array);
    X1 = [X, fliplr(X)];
    color_array = {'r','b','k','c','m','g'};%equal to at least the number of coil models
    marker_array = {'*','o','s','d','p','h'};%equal to at least the number of coil models
    for cx=1:length(coil_models)
        coil_name = coil_models{cx};
        if strcmpi(data_type,'v')
            inBetween = [min(Errors.min_error_GVE{cx}*100,[],1), max(fliplr(Errors.max_error_GVE{cx}*100),[],1)];
            fill(current_axis,X1, inBetween, color_array{cx},'HandleVisibility','off','FaceAlpha',0.2);hold on
            plot(current_axis,Errors.mean_GVE{cx}*100,['-',marker_array{cx},color_array{cx}],'Linewidth',2,'DisplayName',['Mean GVE - ',coil_name]);
            ylabel(current_axis,'GVE (%)')
        elseif strcmpi(data_type,'m')
            inBetween = [min(Errors.min_error_GME{cx}*100,[],1), max(fliplr(Errors.max_error_GME{cx}*100),[],1)];
            fill(current_axis,X1, inBetween, color_array{cx},'HandleVisibility','off','FaceAlpha',0.2);hold on
            plot(current_axis,Errors.mean_GME{cx}*100,['-',marker_array{cx},color_array{cx}],'Linewidth',2,'DisplayName',['Mean GME - ',coil_name]);
            ylabel(current_axis,'GME (%)')
        end
    end
end

function [] = plot_box(current_axis,Errors,coil_models)
    color_array = {'r','b','k','c','m','g'};%equal to at least the number of coil models
    marker_array = {'*','o','s','d','p','h'};%equal to at least the number of coil models
    for cx=1:length(coil_models)
        coil_name = coil_models{cx};
        if strcmpi(data_type,'v')
            boxplot(current_axis,Errors.error_GVE'*100,Errors.Label_GVE,'colors','b'); hold on
            plot(current_axis,Errors.mean_GVE*100,['-',marker_array{cx},color_array{cx}],'Linewidth',2,'DisplayName',['Mean GVE - ',coil_name]);
            ylabel(current_axis,'GVE (%)')
        elseif strcmpi(data_type,'m')
            boxplot(current_axis,Errors.error_GME'*100,Errors.Label_GME,'colors','r');hold on
            plot(current_axis,Errors.mean_GME*100,['-',marker_array{cx},color_array{cx}],'Linewidth',2,'DisplayName',['Mean GME - ',coil_name]);
            ylabel(current_axis,'GME (%)')
        end
    end
end

function [] = plot_range(current_axis,Mode_Array,Errors,data_type,coil_models)
    color_array = {'r','b','k','c','m','g'};%equal to at least the number of coil models
    marker_array = {'*','o','s','d','p','h'};%equal to at least the number of coil models
    for cx=1:length(coil_models)
        coil_name = coil_models{cx};
        if strcmpi(data_type,'v') 
            errorbar(current_axis,1:length(Mode_Array),mean(Errors.mean_GVE{cx}*100,1),mean(Errors.mean_GVE{cx}*100,1)-min(Errors.min_error_GVE{cx},[],1)*100,max(Errors.max_error_GVE{cx}(9:16,:),[],1)*100-mean(Errors.mean_GVE{cx}*100,1),['-',marker_array{cx},color_array{cx}],'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',color_array{cx},'LineWidth',2,'DisplayName',coil_name); hold on
            ylabel(current_axis,'GVE (%)');
        elseif strcmpi(data_type,'m')
            errorbar(current_axis,1:length(Mode_Array),mean(Errors.mean_GME{cx}*100,1),mean(Errors.mean_GME{cx}*100,1)-min(Errors.min_error_GME{cx},[],1)*100,max(Errors.max_error_GME{cx}(9:16,:),[],1)*100-mean(Errors.mean_GME{cx}*100,1),['-',marker_array{cx},color_array{cx}],'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',color_array{cx},'LineWidth',2,'DisplayName',coil_name); hold on
            ylabel(current_axis,'GME (%)')
        end
    end
end

function [] = plot_mean(current_axis,Mode_Array,Errors,data_type,coil_models)
    line_style = {'-','--','-.'};
    color_array = {'r','b','k','c','m','g'};%equal to at least the number of coil models
    marker_array = {'*','o','s','d','p','h'};%equal to at least the number of coil models
    for cx=1:length(coil_models)
        coil_name = coil_models{cx};
        if strcmpi(data_type,'v')
            plot(current_axis,1:length(Mode_Array),mean(Errors.mean_GVE{cx}*100,1),[line_style{cx},marker_array{cx},color_array{cx}],'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color_array{cx},'LineWidth',2,'DisplayName',coil_name); hold on
            ylabel(current_axis,'mean GVE (%)');
        elseif strcmpi(data_type,'m')
            plot(current_axis,1:length(Mode_Array),mean(Errors.mean_GME{cx}*100,1),[line_style{cx},marker_array{cx},color_array{cx}],'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color_array{cx},'LineWidth',2,'DisplayName',coil_name); hold on
            ylabel(current_axis,'mean GME (%)');
        end
    end
end