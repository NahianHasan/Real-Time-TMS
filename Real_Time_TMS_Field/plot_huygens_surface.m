function [current_axis] = plot_huygens_surface(current_axis,msh_file,msh_file_read_fcn,mapping_surface,m2m_dir,...
                                        output_folder,FEMORD,NModes,subject_folder,show_huygens_points,...
                                                                        view_angles,export_fig,save_file)
    hold(current_axis,'on')
    M = matfile(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],...
                                                [subject_folder,'_FEM_',num2str(FEMORD),'.mat']));
    p = M.p;te2p = M.te2p;teid = M.teid;rs = M.rs;surf_point_ind = M.surf_point_ind;
    pp = p;pp(:,surf_point_ind) = M.extruded_scalp_points';
    [GM_tri,~]=surftri(p',te2p(:,teid)');
    [scalp_tri,~]=surftri(p',te2p');
    if exist(fullfile(m2m_dir,'mri2mesh_log.html'),'file')
        model_creation_tool = 'mri2mesh';
    elseif exist(fullfile(m2m_dir,'headreco_log.html'),'file')
        model_creation_tool = 'headreco';
    end
    [~,M1] = get_mapping_surface(msh_file,msh_file_read_fcn,mapping_surface,m2m_dir,model_creation_tool);
    trisurf(scalp_tri,pp(1,:)',pp(2,:)',pp(3,:)','edgealpha',0,'facealpha',0.5,'facecolor',[150 114 100]/255,'DisplayName',"Huygens' Surface",'Parent',current_axis); 
    trisurf(M1.triangles,M1.nodes(:,1),M1.nodes(:,2),M1.nodes(:,3),'edgealpha',0,'facealpha',1,'facecolor','r','Parent',current_axis,'DisplayName','ROI');
    if show_huygens_points
        prc = 0.2;
        random_points = randperm(size(rs,2),ceil(size(rs,2)*1));
        scatter3(current_axis,rs(1,random_points),rs(2,random_points),rs(3,random_points),0.2,'b','filled','DisplayName',...
                                                                                                'Huygens Surface Nodes');
    end
    axis(current_axis,'equal');
    light(current_axis)
    lighting(current_axis,'gouraud');
    c = camlight(current_axis,'headlight');      % Create light
    set(c,'style','infinite');          % Set style
    h = rotate3d;                 % Create rotate3d-handle
    h.ActionPostCallback = @RotationCallback; % assign callback-function
    h.Enable = 'on';              % no need to click the UI-button
    grid(current_axis,'off')
    LGD = legend(current_axis);
    LGD.Location = 'northoutside';
    xlabel(current_axis,'x');ylabel(current_axis,'y');zlabel(current_axis,'z');
    axis(current_axis,'off')
    set(gca,'FontSize',12,'fontweight','bold','linewidth',2)
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
            'ContentType','vector','Resolution',1200,'BackgroundColor','w')
        saveas(current_axis,fullfile([save_file,'.fig']))
    end
    % Sub function for callback
    function RotationCallback(~,~)
        c = camlight(c,'headlight');
    end
end