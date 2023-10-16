function [current_axis] = plot_electric_current_modes(current_axis,msh_file,msh_file_read_fcn,mapping_surface,m2m_dir,...
                                                FEMORD,NModes,output_folder,subject_folder,...
                                                Mode_ID,view_angles,export_fig,randomized,save_file)
	current_mode = load(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],...
                            ['Ax_',num2str(Mode_ID),'.mat'])).Ai;
    M = matfile(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],...
                                                [subject_folder,'_FEM_',num2str(FEMORD),'.mat']));
    p = M.p;te2p = M.te2p;rs = M.rs;
    [scalp_tri,~]=surftri(p',te2p'); 
    if exist(fullfile(m2m_dir,'mri2mesh_log.html'),'file')
        model_creation_tool = 'mri2mesh';
    elseif exist(fullfile(m2m_dir,'headreco_log.html'),'file')
        model_creation_tool = 'headreco';
    end
    [~,M1] = get_mapping_surface(msh_file,msh_file_read_fcn,mapping_surface,m2m_dir,model_creation_tool);
    current_mode = reshape(current_mode,[7,numel(rs)/3,1]);
    electrice_currents = current_mode(1:3,:);
    if randomized
        prc = 0.2;
        rand_id = randperm(size(electrice_currents,2),ceil(prc*size(electrice_currents,2)));
        rs = rs(:,rand_id);
        electrice_currents = electrice_currents(:,rand_id);
    end
    %magnetic_currents = current_mode(4:6,:);
    hold(current_axis,'on');
    trisurf(scalp_tri,p(1,:)',p(2,:)',p(3,:)','edgealpha',0,'facealpha',0.2,'facecolor',[150 114 100]/255,'Parent',current_axis);
    trisurf(M1.triangles,M1.nodes(:,1),M1.nodes(:,2),M1.nodes(:,3),'edgealpha',0,'facealpha',1,'facecolor',[150 114 100]/255,'Parent',current_axis);
    quiver3(current_axis,rs(1,:),rs(2,:),rs(3,:),electrice_currents(1,:),electrice_currents(2,:),electrice_currents(3,:),5,'r','AutoScale','on','AutoScaleFactor',2)
    %quiver3(current_axis,rs(1,:),rs(2,:),rs(3,:),magnetic_currents(1,:),magnetic_currents(2,:),magnetic_currents(3,:),'b')
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
        exportgraphics(current_axis,fullfile([save_file,'_',num2str(Mode_ID),'.png']),...
            'ContentType','vector','Resolution',1200,'BackgroundColor','w')
        saveas(current_axis,fullfile([save_file,'_',num2str(Mode_ID),'.fig']))
    end
    % Sub function for callback
    function RotationCallback(~,~)
        c = camlight(c,'headlight');
    end
end