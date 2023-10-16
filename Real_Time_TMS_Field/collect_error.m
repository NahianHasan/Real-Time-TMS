function [Errors] = collect_error(grid_spacing,Mode_Array,FEMORD,output_folders,hardware,coil_models)
    percentile = 100;
    error_GVE = {};
    mean_error_GVE = {};
    error_GME = {};
    mean_error_GME = {};
    min_error_GVE = {}; max_error_GVE = {};
    min_error_GME = {}; max_error_GME = {};
    Label_GVE = {};Label_GME = {};
    for cx = 1:length(coil_models)
        error_GVE{cx} = [];
        mean_error_GVE{cx} = [];
        error_GME{cx} = [];
        mean_error_GME{cx} = [];
        min_error_GVE{cx} = []; max_error_GVE{cx} = [];
        min_error_GME{cx} = []; max_error_GME{cx} = [];
        Label_GVE{cx} = [];Label_GME{cx} = [];
    end
    for cx=1:length(coil_models)
        coil_model = coil_models{cx};
        for jx=1:length(output_folders)
            output_folder = output_folders{jx};
            for ix=1:length(Mode_Array)
                NModes = Mode_Array(ix);
                if strcmp(hardware,'gpu')
                    W = load(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Error_grid_',num2str(grid_spacing),'_Modes_',num2str(NModes),'_coil_',coil_model,'_gpu.mat']));
                elseif strcmp(hardware,'cpu')
                    W = load(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Error_grid_',num2str(grid_spacing),'_Modes_',num2str(NModes),'_coil_',coil_model,'_cpu.mat']));
                end
                %load GVEs
                temp_error = W.Errors_GVE(:,:,:);
                temp_prc = prctile(temp_error,percentile,1);
                for kx=1:length(temp_prc)
                    temp = temp_error(temp_error(:,kx)<=temp_prc(kx),kx);
                    error_GVE{cx} = [error_GVE{cx};temp];
                    min_error_GVE{cx} = [min_error_GVE{cx},min(temp)];
                    max_error_GVE{cx} = [max_error_GVE{cx},max(temp)];
                    Label_GVE{cx} = [Label_GVE{cx}, repmat({num2str(NModes)},1,length(temp))];
                    mean_error_GVE{cx} = [mean_error_GVE{cx}, mean(temp)];
                end
                %load GMEs
                temp_error1 = W.Errors_GVE1(:,:,:);
                temp_prc1 = prctile(temp_error1,percentile,1);
                for kx=1:length(temp_prc1)
                    temp = temp_error1(temp_error1(:,kx)<=temp_prc1(kx),kx);
                    error_GME{cx} = [error_GME{cx};temp];
                    min_error_GME{cx} = [min_error_GME{cx},min(temp)];
                    max_error_GME{cx} = [max_error_GME{cx},max(temp)];
                    Label_GME{cx} = [Label_GME{cx}, repmat({num2str(NModes)},1,length(temp))];
                    mean_error_GME{cx} = [mean_error_GME{cx}, mean(temp)];
                end
            end
        end
        mean_error_GVE{cx} = reshape(mean_error_GVE{cx},length(Mode_Array),[])';
        mean_error_GME{cx} = reshape(mean_error_GME{cx},length(Mode_Array),[])';
        min_error_GVE{cx} = reshape(min_error_GVE{cx},length(Mode_Array),[])';
        min_error_GME{cx} = reshape(min_error_GME{cx},length(Mode_Array),[])';
        max_error_GVE{cx} = reshape(max_error_GVE{cx},length(Mode_Array),[])';
        max_error_GME{cx} = reshape(max_error_GME{cx},length(Mode_Array),[])';
    end
    
    %create a structure to be returned
    Errors = struct();
    Errors.error_GVE = error_GVE;
    Errors.error_GME = error_GME;
    Errors.mean_error_GVE = mean_error_GVE;
    Errors.mean_error_GME = mean_error_GME;
    Errors.min_error_GVE = min_error_GVE;
    Errors.max_error_GVE = max_error_GVE;
    Errors.min_error_GME = min_error_GME;
    Errors.max_error_GME = max_error_GME;
    Errors.Label_GVE = Label_GVE;
    Errors.Label_GME = Label_GME;
end