function [current_axis] = plot_coil_model(current_axis,real_time_code_path,msh_file,msh_file_read_fcn,m2m_dir,mapping_surface,...
                                            coil_model_file,Transformation,view_angles,export_fig,save_file)
    if exist(fullfile(m2m_dir,'mri2mesh_log.html'),'file')
        model_creation_tool = 'mri2mesh';
    elseif exist(fullfile(m2m_dir,'headreco_log.html'),'file')
        model_creation_tool = 'headreco';
    end
    [~,M] = get_mapping_surface(msh_file,msh_file_read_fcn,mapping_surface,m2m_dir,model_creation_tool);
    [rcoil,~,~,~,coil_tri,~] = load_coil_model(real_time_code_path,coil_model_file,0);
    rcoil2 = rcoil';
    rcoil2(4,:) = 1; rcoil2 = Transformation*rcoil2;
    hold(current_axis,'on');
    trisurf(M.scalp_tri,M.original_msh_nodes(:,1),M.original_msh_nodes(:,2),M.original_msh_nodes(:,3),'edgealpha',0,'facealpha',0.1,'facecolor',[184 115 51]/256,'HandleVisibility','off','DisplayName','','Parent',current_axis);
    trisurf(M.triangles,M.nodes(:,1),M.nodes(:,2),M.nodes(:,3),'edgealpha',0,'facealpha',1,'facecolor',[150 114 100]/255,'HandleVisibility','off','DisplayName','','Parent',current_axis);
    trisurf(coil_tri,rcoil2(1,:)',rcoil2(2,:)',rcoil2(3,:)','edgealpha',0.1,'facealpha',1,'facecolor',[183,119,41]/255,'HandleVisibility','off','DisplayName','','Parent',current_axis);
    %scatter3(rcoil2(1,:),rcoil2(2,:),rcoil2(3,:),10,[150 114 100]/255,'filled','HandleVisibility','off','DisplayName','');hold off;
    light(current_axis)
    lighting(current_axis,'gouraud');
    c = camlight(current_axis,'headlight');      % Create light
    set(c,'style','infinite');          % Set style
    h = rotate3d;                 % Create rotate3d-handle
    h.ActionPostCallback = @RotationCallback; % assign callback-function
    h.Enable = 'on';              % no need to click the UI-button
    grid(current_axis,'off');
    axis(current_axis,'equal')
    axis(current_axis,'off')
    xlabel(current_axis,'X');ylabel(current_axis,'Y');zlabel(current_axis,'Z');
    set(current_axis,'FontSize',12,'fontweight','bold','linewidth',2)
    if view_angles(1)
        view_angle = 2;
    end
    if view_angles(2)
        view_angle = [90 0];
    end
    if view_angles(3)
        view_angle = [0 90];
    end
    if view_angles(4)
        view_angle = 3;
    end
    view(current_axis,view_angle);
    if export_fig
        pathparts = strsplit(save_file,filesep);
        save_dir = fullfile('/',pathparts{1:end-1});
        if ~exist(save_dir,'dir')
            mkdir(save_dir);
        end
        exportgraphics(current_axis,fullfile([save_file,'.png']),...
            'ContentType','vector','Resolution',300,'BackgroundColor','w')
        saveas(current_axis,fullfile([save_file,'.fig']))
    end
    
    % Sub function for callback
    function RotationCallback(~,~)
        c = camlight(c,'headlight');
    end
end