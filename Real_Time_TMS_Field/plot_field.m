function [current_axis] = plot_field(current_axis,msh_file,msh_file_read_fcn,mapping_surface,m2m_dir,...
                                                E_field,view_angles,export_fig,save_file)
    if exist(fullfile(m2m_dir,'mri2mesh_log.html'),'file')
        model_creation_tool = 'mri2mesh';
    elseif exist(fullfile(m2m_dir,'headreco_log.html'),'file')
        model_creation_tool = 'headreco';
    end
    [~,M] = get_mapping_surface(msh_file,msh_file_read_fcn,mapping_surface,m2m_dir,model_creation_tool);
    E_field = reshape(E_field,3,[]);
    E_field_norm = vecnorm(E_field,2,1);
    E_field_norm_mapped = pdeintrp(M.nodes',M.triangles',E_field_norm');
    trisurf(M.triangles,M.nodes(:,1),M.nodes(:,2),M.nodes(:,3),E_field_norm_mapped,'edgealpha',0,'facealpha',1,'Parent',current_axis);
    CX = colorbar(current_axis);
    CX.Limits = [min(E_field_norm_mapped) max(E_field_norm_mapped)];
    CX.Ticks = [min(E_field_norm_mapped):(max(E_field_norm_mapped)-min(E_field_norm_mapped))/5:max(E_field_norm_mapped)];
    CX.TickLabels = compose('%9.2f',CX.Ticks); 
    CX.Title.String = 'V/m';
    light(current_axis)
    lighting(current_axis,'gouraud');
    c = camlight(current_axis,'headlight');      % Create light
    set(c,'style','infinite');          % Set style
    h = rotate3d;                 % Create rotate3d-handle
    h.ActionPostCallback = @RotationCallback; % assign callback-function
    h.Enable = 'on';              % no need to click the UI-button
    grid(current_axis,'off');
    xlabel(current_axis,'X');ylabel(current_axis,'Y');zlabel(current_axis,'Z');
    axis(current_axis,'off')
    axis(current_axis,'equal')
    set(current_axis,'FontSize',20,'fontweight','bold','linewidth',2)
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