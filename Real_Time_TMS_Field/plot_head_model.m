function [current_axis] = plot_head_model(current_axis,msh_file,msh_file_read_fcn,Tissue,view_angles,export_fig,save_file)
    %Tissue = [WM,GM,CSF,skull,scalp]; order from inside to outside of head
    Tissue_Labels = {'WM','GM','CSF','Skull','Scalp'};
    [p,te2p,~,reg,~] = load_msh_data(msh_file,msh_file_read_fcn);
    outermost_tissue_reg_id = max(find(Tissue==1));
    [tri,~]=surftri(p',te2p(:,reg==outermost_tissue_reg_id)'); 
    hold(current_axis,'on');
    trisurf(tri,p(1,:)',p(2,:)',p(3,:)','edgealpha',0,'facealpha',1,'facecolor',[150 114 100]/255,'DisplayName','','Parent',current_axis);        
    light(current_axis)
    lighting(current_axis,'gouraud');
    c = camlight(current_axis,'headlight');      % Create light
    set(c,'style','infinite');          % Set style
    h = rotate3d;                 % Create rotate3d-handle
    h.ActionPostCallback = @RotationCallback; % assign callback-function
    h.Enable = 'on';              % no need to click the UI-button
    grid(current_axis,'off');
    box(current_axis,'on')
    axis(current_axis,'equal')
    xlabel(current_axis,'X');ylabel(current_axis,'Y');zlabel(current_axis,'Z');
    LGD = legend(current_axis,Tissue_Labels{outermost_tissue_reg_id});
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