function [] = ground_truth_job(real_time_code_path,msh_file,msh_file_read_fcn,msh_file_read_fcn_location,m2m_dir,coil_model_file,FEMORD,mapping_surface,output_folder,Transformations,run_mode,cluster_parameter_file,num_simulations)
    if strcmpi(run_mode,'parallel')
        [cluster_params] = collect_cluster_parameters(cluster_parameter_file);
        matlab_module_version = cluster_params{10};
        gt_cluster_name = cluster_params{11};
        gt_cpu = cluster_params{12};
        gt_max_walltime = cluster_params{13};
        
        t = strsplit(coil_model_file,filesep);coil_model = t{end}(1:end-4);
        data_path = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Random_Coil_Placement_IDs_',coil_model,'.mat']);
        save(data_path,'Transformations','-v7.3');
        if ~exist(fullfile('.','slurm_output'),'dir')
            mkdir(fullfile('.','slurm_output'))
        end
        parallel_bash_script_name = fullfile(real_time_code_path,'ground_truth_parallel_run_setup.sh');
        if ispc
            conversion_command = 'unix2dos';
        else
            conversion_command = 'dos2unix';
        end
        cmdStr = [conversion_command ' ' parallel_bash_script_name];
        system(cmdStr);
        cmdStr = [conversion_command ' ' fullfile(real_time_code_path,'ground_truth_cluster_run.sh')];
        system(cmdStr);
        cmdStr = ['chmod +x ' fullfile(real_time_code_path,'ground_truth_cluster_run.sh')];
        system(cmdStr);
        cmdStr = ['bash' ' ' parallel_bash_script_name ' ' real_time_code_path ' ' msh_file ' '...
                    msh_file_read_fcn ' ' msh_file_read_fcn_location ' ' m2m_dir ' ' coil_model_file ' ' num2str(FEMORD) ' '...
                    mapping_surface ' ' output_folder ' ' num2str(num_simulations) ' ' gt_cluster_name ' ' gt_cpu ' '...
                    gt_max_walltime ' ' matlab_module_version];
        system(cmdStr);
        disp([newline,'The results will be stored inside the output folder',newline])
    elseif strcmpi(run_mode,'SERIAL')
        for ix = 1:num_simulations
            mode = ix;%a random identifier for the simulation;
            Transformation = Transformations(:,:,ix);
            [~,~,~] = ground_truth(real_time_code_path,mode,msh_file,msh_file_read_fcn,msh_file_read_fcn_location,m2m_dir,coil_model_file,FEMORD,mapping_surface,output_folder,Transformation);
        end
        disp([newline,'The results are stored inside the output folder',newline]);
    else
        disp("Please specify whether to run the offline stage in 'parallel' or 'serial'")
    end
end