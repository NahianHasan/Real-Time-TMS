function [output_folder] = compile(real_time_code_path,msh_file_read_fcn_location,output_directory,subject_folder)
    addpath(real_time_code_path)
    addpath(fullfile(real_time_code_path,'matlab'));
    addpath(msh_file_read_fcn_location);
    output_folder=fullfile(output_directory,subject_folder);
    if ~exist(output_folder,'dir')
        mkdir(output_folder);
    end
end