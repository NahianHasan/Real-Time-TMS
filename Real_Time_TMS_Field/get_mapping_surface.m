function [GM_teid,GM_msh] = get_mapping_surface(msh_file,msh_file_read_fcn,mapping_surface,m2m_dir,model_creation_tool)
    [p,te2p,conductivity,~,GM_msh] = load_msh_data(msh_file,msh_file_read_fcn);
    if strcmpi(mapping_surface, 'WM')% WM surface
        xx = [];
        xx(conductivity==.1260)=1; % White
        [tri,GM_teid]=surftri(p',te2p(:,xx(:)==1)'); %WM surface
        GM_msh.triangles = tri;
        GM_msh.triangle_regions = ones(size(tri,1),1).*1;
    elseif strcmpi(mapping_surface, 'GM')% GM+WM surface
        xx = [];
        xx(conductivity==.1260)=1; % White
        xx(conductivity==.2750)=1; % Grey
        [tri,GM_teid]=surftri(p',te2p(:,xx(:)==1)'); % GM/WM surface
        GM_msh.triangles = tri;
        GM_msh.triangle_regions = ones(size(tri,1),1).*2;
    elseif strcmpi(mapping_surface, 'GMM')% GM middle surface
        teid=1:numel(te2p)/4;
        teid=teid(conductivity==0.2750 | conductivity==0.1260);
        [~,GMM_msh] = load_GM_mid_Layer(m2m_dir,model_creation_tool);
        GM_msh.nodes = GMM_msh.nodes./1000;
        GM_msh.triangles = GMM_msh.triangles;
        GM_msh.triangle_regions = GMM_msh.triangle_regions;
        v1 = GM_msh.nodes(GM_msh.triangles(:,2),:)-GM_msh.nodes(GM_msh.triangles(:,1),:);
        v2 = GM_msh.nodes(GM_msh.triangles(:,3),:)-GM_msh.nodes(GM_msh.triangles(:,1),:);
        SA = sum(0.5*vecnorm(cross(v1,v2),2,2),'all')/1E6;
        warning('off');TR = triangulation(te2p(:,teid)',p');warning('on')
        GM_teid = pointLocation(TR,GM_msh.nodes);
        GM_teid(isnan(GM_teid)) = 1;
        clear TR
    end
    GM_msh.tetrahedra = [];
    GM_msh.tetrahedron_regions = [];
    GM_msh.conductivity = [];
end

