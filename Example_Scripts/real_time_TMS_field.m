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
msh_file='/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.2-Realtime_TMS/Code_Github/Data/ernie_mri2mesh/ernie.msh';
msh_file_read_fcn = 'mesh_load_gmsh4';%useful for '.msh' msh-files. Also make sure the function is in the matlab search path. For '.mat' msh-files, it's not necessary.
msh_file_read_fcn_location = fullfile('/home/wang3007-1/SimNIBS-3.2/matlab');
m2m_dir = '/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.2-Realtime_TMS/Code_Github/Data/ernie_mri2mesh/m2m_ernie';
FEMORD=1;%FEM order = 1,2, or, 3
mapping_surface = 'GMM';%options: 'WM','GM'=cortex, 'GMM'=middle grey matter surface
grid_spacing=0.004;%interpolation grid spacing. Recommended size is <=0.004 (i.e., <= 4mm). Used for primary field calculation.
output_directory = '/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.2-Realtime_TMS/Code_Github/Solutions';%output directory for the results
subject_folder = 'ernie_mri2mesh';%name of the subject
coil_model_file = '/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.2-Realtime_TMS/Code_Github/Coil_Models/coil_fig-8.mat';
run_mode='serial';%options = 'serial','parallel' (for HPC clusters);
%If parallel, provide the cluster parameters in a separate csv file (cluster_parameters.csv) (compatible with slurm scripting)
cluster_parameter_file = '/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.2-Realtime_TMS/Code_Github/Example_Scripts/cluster_parameters.csv';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TMS E-field Prediction for Multiple Coil Placements for Analysis
[output_folder] = compile(real_time_code_path,msh_file_read_fcn_location,...
                            output_directory,subject_folder);
%set_up the real_time stage
[Q,Ax,params] = real_time_stage_setup(msh_file,msh_file_read_fcn,real_time_code_path,...
        NModes,m2m_dir,FEMORD,grid_spacing,mapping_surface,output_folder,...
        subject_folder,coil_model_file);
%provide a 4*4 transformation matrix for the coil placement relative to the
%head. In real-time scenario, this transformation matrix comes from the
%neuronavigation system.
Transformation = rand([4,4]);
%predict the E-field in real time. For faster calculation. If any gpu is
%available, the returned Efield is a gpuArray (still inside the GPU) for
%fast rendering to the real-time neuro navigation. The field will be 
%accessible until the gpu is reset. Note: the execution time
%also includes the communication time of transformation matrix and its
%inverse from host to device.
[Efield] = real_time_field_calculation_single_placements(Q,Ax,params,Transformation);
