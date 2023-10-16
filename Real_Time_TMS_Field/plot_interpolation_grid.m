function [current_axis] = plot_interpolation_grid(current_axis,real_time_code_path,msh_file,msh_file_read_fcn,...
                           mapping_surface,m2m_dir,output_folder,FEMORD,NModes,subject_folder,...
                           coil_model_file,view_angles,export_fig,save_file,Transformation)
    [rcoil,kcoil,~,~,coil_tri,~] = load_coil_model(real_time_code_path,coil_model_file,0);
    
    M = matfile(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],...
                                                [subject_folder,'_FEM_',num2str(FEMORD),'.mat']));
    p = M.p;te2p = M.te2p;
    surf_point_ind = M.surf_point_ind;
    pp = p;pp(:,surf_point_ind) = M.extruded_scalp_points';
    [scalp_tri,~]=surftri(p',te2p');
    pp(4,:) = 1;
    pp_rot = Transformation\pp;
    if exist(fullfile(m2m_dir,'mri2mesh_log.html'),'file')
        model_creation_tool = 'mri2mesh';
    elseif exist(fullfile(m2m_dir,'headreco_log.html'),'file')
        model_creation_tool = 'headreco';
    end
    [~,M1] = get_mapping_surface(msh_file,msh_file_read_fcn,mapping_surface,m2m_dir,model_creation_tool);
    nodes = M1.nodes';
    nodes(4,:) = 1;
    nodes_rot = Transformation\nodes;
    %generating a coarse grid of 10 mm just for visualization, actual grid
    %size should be <= 4mm as recommended in the paper.
    [~,~,~,~,X,Y,Z,~,~,~] = volumetric_grid(p,0.01,rcoil,kcoil,1);
    hold(current_axis,"on");
    trisurf(scalp_tri,pp_rot(1,:)',pp_rot(2,:)',pp_rot(3,:)','edgealpha',0,'facealpha',0.5,'facecolor',[150 114 100]/255,'DisplayName',"Huygens' Surface",'Parent',current_axis); 
    trisurf(M1.triangles,nodes_rot(1,:)',nodes_rot(2,:)',nodes_rot(3,:)','edgealpha',0,'facealpha',1,'facecolor','r','Parent',current_axis,'DisplayName','ROI');
    trisurf(coil_tri,rcoil(:,1),rcoil(:,2),rcoil(:,3),'edgealpha',0,'facealpha',1,'facecolor',[183,119,41]/255,'DisplayName','TMS Coil','Parent',current_axis);    
    scatter3(current_axis,reshape(X,1,[]),reshape(Y,1,[]),reshape(Z,1,[]),2,'g','filled','DisplayName','Grid Locations');
    axis(current_axis,'equal');
    light(current_axis)
    lighting(current_axis,'gouraud');
    c = camlight(current_axis,'headlight');      % Create light
    set(c,'style','infinite');          % Set style
    h = rotate3d;                 % Create rotate3d-handle
    h.ActionPostCallback = @RotationCallback; % assign callback-function
    h.Enable = 'on';              % no need to click the UI-button
    grid(current_axis,'off')
    xlabel(current_axis,'x');ylabel(current_axis,'y');zlabel(current_axis,'z');
    axis(current_axis,'off')
    set(gca,'FontSize',12,'fontweight','bold','linewidth',2)
    if view_angles(1)
        view_angle = 2;
    end
    if view_angles(2)
        view_angle = [180 0];
    end
    if view_angles(3)
        view_angle = [0 -90];
    end
    if view_angles(4)
        view_angle = 3;
    end
    view(current_axis,view_angle);
    LGD = legend(current_axis);
    LGD.Location = 'northoutside';
    if export_fig
        pathparts = strsplit(save_file,filesep);
        save_dir = fullfile('/',pathparts{1:end-1});
        if ~exist(save_dir,'dir')
            mkdir(save_dir);
        end
        exportgraphics(current_axis,fullfile([save_file,'.png']),...
            'ContentType','vector','Resolution',1200,'BackgroundColor','w')
        saveas(current_axis,fullfile([save_file,'.fig']))
    end
    % Sub function for callback
    function RotationCallback(~,~)
        c = camlight(c,'headlight');
    end
end