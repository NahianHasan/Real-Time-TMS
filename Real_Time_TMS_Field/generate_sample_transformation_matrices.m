function [Anor,Transformation,coil_angles,coil_positions,coil_position_ids,combinations] = generate_sample_transformation_matrices(msh_file,m2m_dir,coil_model_file,num_simulations,simnibs_installation_dir,th_hair)
    addpath(simnibs_installation_dir);
    addpath(fullfile(simnibs_installation_dir,'matlab'));
    warning('off');[pp,Anor,~] = generate_sample_coil_placement(msh_file,m2m_dir,th_hair);warning('on');
    combinations = randperm(size(Anor,3)*360,num_simulations);

    Transformation = zeros(4,4,length(combinations));
    coil_angles = zeros(1,length(combinations));
    coil_positions = zeros(3,length(combinations));
    coil_position_ids = zeros(1,length(combinations));
    for kx=1:length(combinations)
        coil_ang = mod(combinations(kx),360);
        coil_pos = floor(combinations(kx)/360);
        if coil_ang == 0
            coil_ang = 360;
        else
            coil_pos = coil_pos + 1;
        end
        rot_mat = [cos(coil_ang/180*pi), sin(coil_ang/180*pi), 0, 0;-sin(coil_ang/180*pi), cos(coil_ang/180*pi), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
        Transformation(:,:,kx) = Anor(:,:,coil_pos)*rot_mat;
        coil_angles(kx) = coil_ang;
        coil_positions(:,kx) = pp(coil_pos,:)';
        coil_position_ids(kx) = coil_pos;
    end
end