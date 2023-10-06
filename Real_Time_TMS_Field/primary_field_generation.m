function [Eout,Hout,rcoil,kcoil,X,Y,Z] = primary_field_generation(msh_file,msh_file_read_fcn,real_time_code_path,coil_model_file,grid_spacing,FEMORD,output_folder)
    [p,~,~,~] = load_msh_data(msh_file,msh_file_read_fcn);
    [rcoil,kcoil,coil_model,~,~,~] = load_coil_model(real_time_code_path,coil_model_file,1);
    [Eout,Hout,rcoil,kcoil,X,Y,Z,coil_thickness,min_coil,del_z] = volumetric_grid(p,grid_spacing,rcoil,kcoil);
    save(fullfile(output_folder,['FEM_',num2str(FEMORD)],['grid_fields_',num2str(grid_spacing),'_',coil_model,'.mat']),'Eout','Hout','rcoil','kcoil','coil_thickness','min_coil','X','Y','Z','del_z','-v7.3');
end
