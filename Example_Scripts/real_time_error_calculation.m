clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defining Paths
real_time_code_path = '/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.2-Realtime_TMS/Code_Github/Real_Time_TMS_Field';
addpath(real_time_code_path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter Set-up
NModes=400;%number of modes
msh_file='/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.2-Realtime_TMS/Code_commercial/data/ernie_mri2mesh/ernie.msh';
msh_file_read_fcn = 'mesh_load_gmsh4';%useful for '.msh' msh-files. Also make sure the function is in the matlab search path. For '.mat' msh-files, it's not necessary.
msh_file_read_fcn_location = fullfile('/home/wang3007-1/SimNIBS-3.2/matlab');
m2m_dir = '/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.2-Realtime_TMS/Code_commercial/data/ernie_mri2mesh/m2m_ernie';
FEMORD=1;%FEM order = 1,2, or, 3
mapping_surface = 'GMM';%options: 'GM'-cortex, 'GMM'-middle grey matter surface
grid_spacing=0.004;%interpolation grid spacing. Recommended size is <=0.004 (i.e., <= 4mm). Used for primary field calculation.
output_directory = '/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.2-Realtime_TMS/Code_Github/solutions_sample';%output directory for the results
subject_folder = 'ernie_mri2mesh';%name of the subject
coil_model_file = '/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.2-Realtime_TMS/Code_Github/Coil_Models/coil_fig-8.mat';
run_mode='serial';%options = 'serial','parallel' (for HPC clusters);
%If parallel, provide the cluster parameters in a separate csv file (cluster_parameters.csv) (compatible with slurm scripting)
cluster_parameter_file = '/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.2-Realtime_TMS/Code_Github/Example_Scripts/cluster_parameters.csv';





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Error Calculation
[output_folder] = compile(real_time_code_path,msh_file_read_fcn_location,output_directory,subject_folder);
%set_up the real_time stage
[Q,Ax,params] = real_time_stage_setup(msh_file,msh_file_read_fcn,...
        real_time_code_path,NModes,m2m_dir,FEMORD,grid_spacing,...
        mapping_surface,output_folder,subject_folder,coil_model_file);
%The ground truth data needs to be calculated beforehand and stored in the
%same output folder.
[Errors_GVE,Errors_GME,Errors_LVE,Errors_LME] = error_calculation(Q,Ax,params,...
                        coil_model_file,FEMORD,NModes,grid_spacing,output_folder);












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
