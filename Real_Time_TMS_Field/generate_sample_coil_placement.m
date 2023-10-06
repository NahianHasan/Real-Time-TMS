function [pp_test,Anor_test,nhatp] = generate_sample_coil_placement(msh_file,m2m_dir,th_hair)
    M = mesh_load_gmsh4(msh_file);
    p = double(M.nodes');
    p = p/1000;%conversion to m scale
    te2p = double(M.tetrahedra');
    scalp_tri=surftri(p',te2p');%select scalp triangular facets

    %Generate some sample coil placements
    use_scalp_projection = 1;
    use_nhat_projection = 0;
    use_msh_smoothing = 1;%Whether to smooth mesh for calculating nhats. 
    %1=smooth EEG-coordinates/projected-scalp-coordinates a bit for nhat calculation
    nhat_normalization = 'angle_mean_normal';
    eeg_mni_source_file = 'EEG10-10_UI_Jurak_2007.csv';    
    sphere_density = 5;
    patch_angle = 60;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(eeg_mni_source_file,'EEG10-10_UI_Jurak_2007.csv')
        eeg_coregistration_exists = 1;
        eeg_csv = [m2m_dir,'/eeg_positions/',eeg_mni_source_file];
        if ~exist(eeg_csv, 'file')
            eeg_coregistration_exists = 0;
            eeg_csv = ['./Real_Time_TMS_Field/ElectrodeCaps_MNI/',eeg_mni_source_file];
        end
    else
        eeg_coregistration_exists = 0;
        eeg_csv = ['./Real_Time_TMS_Field/ElectrodeCaps_MNI/',eeg_mni_source_file];
    end
    F=fopen(eeg_csv,'r');
    tline = fgetl(F);
    EEG_data_points = [];
    EEG_channels = {};
    cnt=1;
    while ischar(tline)
        L=tline;
        t=strsplit(L,',');
        channel = t{end};
        EEG_data_points = [EEG_data_points; [str2double(t{2}), str2double(t{3}), str2double(t{4})]];
        EEG_channels{cnt} = channel;
        cnt = cnt+1;
        tline = fgetl(F);
    end
    fclose(F);
    if eeg_coregistration_exists
        EEG_data_points = EEG_data_points./1000;
        disp('EEG data read')
    else
        %Transform them to subject space
        EEG_data_points = mni2subject_coords(EEG_data_points, m2m_dir); 
        EEG_data_points = EEG_data_points./1000;
        disp('EEG data read')
    end
    
    %generate sample coil positions and corresponding Anor transformations
    [pp_test,Anor_test,~,nhatp]=findSurf_real_head(p,th_hair,EEG_data_points,EEG_channels,scalp_tri,sphere_density,use_scalp_projection,use_nhat_projection,use_msh_smoothing,nhat_normalization,patch_angle);
end