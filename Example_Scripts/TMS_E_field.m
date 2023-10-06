clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defining Paths
real_time_code_path = '/scratch/bell/hasan34/data/Real_Time_TMS/Code_Github/Real_Time_TMS_Field';
addpath(real_time_code_path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter Set-up
NModes=400;%number of modes
msh_file='/scratch/bell/hasan34/data/head_model_data/mri2mesh_models/Updated/BA40/mri2mesh_models/S35/S35.msh';
msh_file_read_fcn = 'mesh_load_gmsh4';%useful for '.msh' msh-files. Also make sure the function is in the matlab search path. For '.mat' msh-files, it's not necessary.
msh_file_read_fcn_location = fullfile('/home/hasan34/SimNIBS-3.2/matlab');
m2m_dir = '/scratch/bell/hasan34/data/head_model_data/mri2mesh_models/Updated/BA40/mri2mesh_models/S35/m2m_S35';
FEMORD=1;%FEM order = 1,2, or, 3
mapping_surface = 'GMM';%options: 'GM'-cortex, 'GMM'-middle grey matter surface
grid_spacing=0.004;%interpolation grid spacing. Recommended size is <=0.004 (i.e., <= 4mm). Used for primary field calculation.
output_directory = '/scratch/bell/hasan34/data/Real_Time_TMS/Code_Github/solutions_sample';%output directory for the results
subject_folder = 'BA40_S35_mri2mesh';%name of the subject
coil_model_file = '/scratch/bell/hasan34/data/Real_Time_TMS/Code_Github/Coil_Models/coil_fig-8.mat';
run_mode='serial';%options = 'serial','parallel' (for HPC clusters);
%If parallel, provide the cluster parameters in a separate csv file (cluster_parameters.csv) (compatible with slurm scripting)
cluster_parameter_file = '/scratch/bell/hasan34/data/Real_Time_TMS/Code_Github/Example_Scripts/cluster_parameters.csv';
simnibs_installation_dir = '/home/hasan34/SimNIBS-3.2';

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TMS E-field Prediction for Multiple Coil Placements for Analysis
[output_folder] = compile(real_time_code_path,msh_file_read_fcn_location,output_directory,subject_folder);
%set_up the real_time stage
[Q,Ax,params] = real_time_stage_setup(msh_file,msh_file_read_fcn,real_time_code_path,NModes,m2m_dir,FEMORD,grid_spacing,mapping_surface,output_folder,subject_folder,coil_model_file);
%generate sample coil placements and transformation matrices based on SimNIBS. You can have your own generator.
num_simulations = 1000;%define the number of placements needed
th_hair = 0.005;%distance from scalp to coil center (in meter).
[~,Transformations,~,~,~,~] = generate_sample_transformation_matrices(msh_file,m2m_dir,coil_model_file,num_simulations,simnibs_installation_dir,th_hair);

%predict the E-field in real time. For faster calculation, set memory_calculation=0;
%The returned fields are in the host (i.e., cpu).
memory_calculation = 0;
[Efields,Time,Memory,data_set_up_time,communication_time,map_teid,map_msh] = real_time_field_calculation_multiple_placements(Q,Ax,params,Transformations,memory_calculation);
