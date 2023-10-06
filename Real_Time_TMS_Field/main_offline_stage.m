function [] = main_offline_stage(real_time_code_path,msh_file_read_fcn_location,NModes,msh_file,msh_file_read_fcn,FEMORD,grid_spacing,output_folder,run_mode,cluster_parameter_file)
    if strcmpi(run_mode,'parallel')
        [cluster_params] = collect_cluster_parameters(cluster_parameter_file);
        cluster_name = cluster_params{1};
        stage_1_cpu = cluster_params{2};
        stage_2_cpu = cluster_params{3};
        stage_3_cpu = cluster_params{4};
        stage_4_cpu = cluster_params{5};
        stage_1_max_walltime = cluster_params{6};
        stage_2_max_walltime = cluster_params{7};
        stage_3_max_walltime = cluster_params{8};
        stage_4_max_walltime = cluster_params{9};
        matlab_module_version = cluster_params{10};
    end
    if ~exist(output_folder,'dir')
        mkdir(output_folder)
    end
    if strcmpi(run_mode,'parallel')
        if ~exist(fullfile('.','slurm_output'),'dir')
            mkdir(fullfile('.','slurm_output'))
        end
        parallel_bash_script_name = fullfile(real_time_code_path,'offline_parallel_run_script.sh');
        if ispc
            conversion_command = 'unix2dos';
        else
            conversion_command = 'dos2unix';
        end
        cmdStr = [conversion_command ' ' parallel_bash_script_name];
        system(cmdStr);
        for stage=1:4
            cmdStr = [conversion_command,' ',fullfile(real_time_code_path,['offline_parallel_stage_',num2str(stage),'.sh'])];
            system(cmdStr);
        end
        cmdStr = ['bash' ' ' parallel_bash_script_name ' ' num2str(NModes) ' ' msh_file ' ' num2str(FEMORD) ' ' ...
                   num2str(grid_spacing) ' ' output_folder ' ' msh_file_read_fcn ' ' stage_1_cpu ' ' ...
                   stage_2_cpu ' ' stage_3_cpu ' ' stage_4_cpu ' ' cluster_name ' ' stage_1_max_walltime ' ' ...
                   stage_2_max_walltime ' ' stage_3_max_walltime ' ' stage_4_max_walltime ' ' ...
                   matlab_module_version ' ' real_time_code_path ' ' msh_file_read_fcn_location];
        system(cmdStr);
        disp([newline,'Modes are being calculated in the cluster',newline])
    elseif strcmpi(run_mode,'serial')
        offline_serial_run_script(msh_file_read_fcn,msh_file_read_fcn_location,NModes,msh_file,FEMORD,output_folder);
    else
        disp("Please specify whether to run the offline stage in 'parallel' or 'serial'")
    end
end
