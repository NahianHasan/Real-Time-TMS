clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defining Paths
real_time_code_path = '/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.2-Realtime_TMS/Code_Github/Real_Time_TMS_Field';
addpath(real_time_code_path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter Set-up
msh_file='/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.2-Realtime_TMS/Code_commercial/data/ernie_mri2mesh/ernie.msh';
msh_file_read_fcn = 'mesh_load_gmsh4';%useful for '.msh' msh-files. Also make sure the function is in the matlab search path. For '.mat' msh-files, it's not necessary.
msh_file_read_fcn_location = fullfile('/home/wang3007-1/SimNIBS-3.2/matlab');
m2m_dir = '/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.2-Realtime_TMS/Code_commercial/data/ernie_mri2mesh/m2m_ernie';
mapping_surface = 'GMM';%options = 'GM','WM','GMM'
output_directory = '/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.2-Realtime_TMS/Code_Github/Solutions';%output directory for the results
coil_model_file = '/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.2-Realtime_TMS/Code_Github/Coil_Models/coil_fig-8.mat';
subject_folder = 'ernie_mri2mesh';%name of the subject
simnibs_installation_dir = '/home/wang3007-1/SimNIBS-3.2';

%%
[output_folder] = compile(real_time_code_path,msh_file_read_fcn_location,output_directory,subject_folder);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Coil Model
% fig1 = figure();
% save_file = fullfile(output_folder,'Results','coil_model');
% [~,Transformations,~,~,~,~] = generate_sample_transformation_matrices(msh_file,m2m_dir,coil_model_file,...
%                                     1000,simnibs_installation_dir,0.005);
% Transformation = Transformations(:,:,1);
% plot_coil_model(gca(fig1),real_time_code_path,msh_file,msh_file_read_fcn,m2m_dir,mapping_surface,...
%                     coil_model_file,Transformation,[1,0,0,0],1,save_file);
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Head Model
% fig2 = figure();
% save_file = fullfile(output_folder,'Results','head_model');            
% plot_head_model(gca(fig2),msh_file,msh_file_read_fcn,[1,0,0,0,0],[1,0,0,0],1,save_file)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Error
% fig3 = figure();
% save_dir = './Results';
% grid_spacing = 0.005;
% FEMORD = 1;
% Mode_Array = [400];%create list of all modes
% output_folders = {output_folder};%create a list of all output folders corresponding to each subject.
% coil_model_files = {coil_model_file};%create a list of all coil model files
% hardware = 'cpu';%define the hardware where the errors were calculated.options='cpu','gpu'
% plot_type = 'mean';%options='box','patch','mean','range'
% data_type = 'v';%options='v' for GVE plot, 'm' for GME plot
% [Errors,Ax] = plot_error(gca(fig3),real_time_code_path,grid_spacing,Mode_Array,FEMORD,hardware,...
%                                         output_folders,save_dir,coil_model_files,plot_type,data_type);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Time
% fig4 = figure();
% save_dir = './Results';
% grid_spacing = 0.005;
% FEMORD = 1;
% Mode_Array = [400];%create list of all modes
% output_folders = {output_folder};%create a list of all output folders corresponding to each subject.
% hardware = 'cpu';%define the hardware where the timings were calculated.options='cpu','gpu'
% plot_type = 'mean';%options='box','patch','mean','range'
% data_type = 'real_time';%options = 'real_time', 'preprocessing_time'
% num_simulations = 1000;%# of simulations.
% [Time,Ax] = plot_time(gca(fig4),real_time_code_path,grid_spacing,Mode_Array,FEMORD,hardware,output_folders,...
%                                                  num_simulations,subject_folders,save_dir,coil_models,plot_type);
                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Memory
% fig5 = figure();
% save_dir = './Results';
% grid_spacing = 0.005;
% FEMORD = 1;
% Mode_Array = [400];%create list of all modes
% output_folders = {output_folder};%create a list of all output folders corresponding to each subject.
% plot_type = 'mean';%options='box','patch','mean','range'
% [Memory,Ax] = plot_memory(gca(fig5),grid_spacing,Mode_Array,FEMORD,output_folders,...
%                                                             save_dir,plot_type,coil_models);
%                                                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot E-field
% fig6 = figure();
% save_file = fullfile(output_folder,'Results','E_field_1');
% %provide the appropriate E-field information based on your simulation
% E_field = Efields(:,1);
% [current_axis] = plot_field(gca(fig6),msh_file,msh_file_read_fcn,mapping_surface,m2m_dir,...
%                                                 E_field,[1,0,0,0],1,save_file);
%                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Huygens Surface
% fig7 = figure();
% save_file = fullfile(output_folder,'Results','Huygens_Surface');      
% [current_axis] = plot_huygens_surface(gca(fig7),msh_file,msh_file_read_fcn,'GMM',m2m_dir,...
%                                         output_folder,FEMORD,NModes,subject_folder,1,...
%                                                                     [1,0,0,0],1,save_file);
%                                     
                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot E-field Modes
% fig8 = figure();
% save_file = fullfile(output_folder,'Results','E_field_mode');
% field_domain = 'GM';%The solution domain of E-field modes. For TMS, 'GM' corresponds to the all elements inside the GM surface.
% [current_axis] = plot_E_field_modes(gca(fig8),output_folder,msh_file,msh_file_read_fcn,field_domain,...
%                                     m2m_dir,FEMORD,NModes,Mode_ID,[1,0,0,0],1,save_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Surface Current Modes
% fig9 = figure();
% save_file = fullfile(output_folder,'Results','Electric_current_mode'); 
% [current_axis] = plot_electric_current_modes(gca(fig9),msh_file,msh_file_read_fcn,'GMM',m2m_dir,...
%                                                 FEMORD,NModes,output_folder,subject_folder,...
%                                                 2,[1,0,0,0],1,1,save_file);
% fig10 = figure();
% save_file = fullfile(output_folder,'Results','Magnetic_current_mode');   
% [current_axis] = plot_magnetic_current_modes(gca(fig10),msh_file,msh_file_read_fcn,'GMM',m2m_dir,...
%                                                 FEMORD,NModes,output_folder,subject_folder,...
%                                                 2,[1,0,0,0],1,1,save_file);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Interpolation Grid
% fig11 = figure();
% save_file = fullfile(output_folder,'Results','Interpolation_Grid'); 
% [~,Transformations,~,~,~,~] = generate_sample_transformation_matrices(msh_file,m2m_dir,coil_model_file,...
%                                     1000,simnibs_installation_dir,0.005);
% Transformation = Transformations(:,:,1);
% [current_axis] = plot_interpolation_grid(gca(fig11),real_time_code_path,msh_file,msh_file_read_fcn,...
%                            mapping_surface,m2m_dir,output_folder,FEMORD,NModes,subject_folder,...
%                            coil_model_file,[0,1,0,0],1,save_file,Transformation);
                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%