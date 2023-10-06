function [Errors_GVE,Errors_GME,Errors_LVE,Errors_LME] = error_calculation(Q,Ax,params,coil_model_file,FEMORD,NModes,grid_spacing,output_folder)
    %collect the ground truth simulation ids
    t = strsplit(coil_model_file,filesep);coil_model = t{end}(1:end-4);
    gt_IDs = [];
    gt_directory = fullfile(output_folder,['FEM_',num2str(FEMORD)],['GT_E_Fields_',coil_model]);
    F = dir(fullfile(gt_directory,'E_org*.mat'));
    for ix=1:length(F)
        t = strsplit(F(ix).name,'_');
        gt_IDs = [gt_IDs,str2double(t{end}(1:end-3))];
    end
    gt_IDs = sort(gt_IDs);
    
    save_file_name = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Error_grid_',num2str(grid_spacing),'_Modes_',num2str(NModes),'_coil_',coil_model,'_gpu.mat']);
    Errors_GVE = [];Errors_GME = [];Errors_LVE = [];Errors_LME = [];
    for jx=1:length(gt_IDs)
        disp(['Loading GT = ',num2str(gt_IDs(jx))])
        try
            M = load(fullfile(gt_directory,['E_org_',num2str(gt_IDs(jx)),'.mat']));
            Efield_gt = M.E_org;
            if size(Efield_gt,1)<3
                Efield_gt = reshape(Efield_gt,3,[]);
            end
            Transformation = M.Transformation;
            Efield_rt = real_time_field_calculation_single_placements(Q,Ax,params,Transformation);
            hardware = 'CPU';
            if numel(params{31})>0
                Efield_rt = gather(Efield_rt);
                hardware = 'GPU';
            end
            
            Efield_rt = reshape(Efield_rt,3,[]);
            Efield_gt_m = vecnorm(Efield_gt,2,1);
            Efield_rt_m = vecnorm(Efield_rt,2,1);
            
            GVE = norm(Efield_rt(:)-Efield_gt(:))/norm(Efield_gt(:));
            GME = norm(Efield_rt_m(:)-Efield_gt_m(:))/norm(Efield_gt_m(:));
            LVE = vecnorm(Efield_rt-Efield_gt,2,1)./max(Efield_gt_m);
            LME = abs(Efield_rt_m-Efield_gt_m)./max(Efield_gt_m);
            
            Errors_GVE = [Errors_GVE;GVE];Errors_GME = [Errors_GME;GME];
            Errors_LVE = [Errors_LVE;LVE];Errors_LME = [Errors_LME;LME];
    
            disp([hardware,' Running - Modes (w) = ',num2str(NModes),'--- ID = ',num2str(jx),...
                        '--- GVE=',num2str(GVE(end)*100),'%, GME=',num2str(GME(end)*100),'%']);            
        catch
            disp('Could not Load');
            continue;
        end
    end
    save(save_file_name,'Errors_GVE','Errors_GME','Errors_LVE','Errors_LME','-v7.3');
    reset(g);
end