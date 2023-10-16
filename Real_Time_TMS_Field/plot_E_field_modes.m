function [current_axis] = plot_E_field_modes(current_axis,output_folder,msh_file,msh_file_read_fcn,field_domain,...
                                            m2m_dir,FEMORD,NModes,Mode_ID,view_angles,export_fig,save_file)
    if exist(fullfile(m2m_dir,'mri2mesh_log.html'),'file')
        model_creation_tool = 'mri2mesh';
    elseif exist(fullfile(m2m_dir,'headreco_log.html'),'file')
        model_creation_tool = 'headreco';
    end
    [mapping_surf_teid,M] = get_mapping_surface(msh_file,msh_file_read_fcn,field_domain,m2m_dir,model_creation_tool);
    field_mode = load(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],...
                            ['Q_',num2str(Mode_ID),'.mat'])).Qi;
    field_mode = reshape(field_mode,3,[]);
    field_mode = field_mode(:,mapping_surf_teid);
    field_mode = vecnorm(field_mode,2,1);
    field_mode = field_mode./max(field_mode);
    trisurf(M.triangles,M.nodes(:,1),M.nodes(:,2),M.nodes(:,3),field_mode,'edgealpha',0,'facealpha',1,'Parent',current_axis);
    CX = colorbar(current_axis);
    CX.Limits = [min(field_mode) max(field_mode)];
    CX.Ticks = [min(field_mode),max(field_mode)];
    %CX.TickLabels = compose('%9.2f',CX.Ticks); 
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