function [GM_mid_centers,m] = load_GM_mid_Layer(m2m_folder,model_creation_tool)
    
    if strcmp(model_creation_tool,'charm')
        % determine what to load
        SIMNIBSDIR = '/home/hasan34/SimNIBS-3.2';
        path_to_avg_surf = fullfile(SIMNIBSDIR, 'resources', 'templates', 'fsaverage_surf');
        path_to_labels = fullfile(SIMNIBSDIR, 'resources', 'templates', 'fsaverage_atlases');
        
        % load subject-specific surfaces
        if exist(fullfile(m2m_folder,'charm_log.html'),'file')
            lh={ fullfile(m2m_folder,'surfaces','lh.central.gii') };
            rh={ fullfile(m2m_folder,'surfaces','rh.central.gii') };
        else
            error(['No .._log.html found in ' m2m_folder '. Unclear whether it was created by charm']);
        end
        for i=1:length(lh)
            if ~exist(lh{i},'file'); error(['could not find ' lh{i}]); end
        end

        for i=1:length(rh)
            if ~exist(rh{i},'file'); error(['could not find ' rh{i}]); end
        end

        % load surfaces
        m=mesh_empty;

        if ~isempty(lh)
            [nodes, triangles] = load_surface(lh{1});
            m.nodes = [m.nodes; nodes];
            m.triangles= [m.triangles; triangles];
            m.triangle_regions= [m.triangle_regions; ones(size(triangles,1),1)];
        end
        if length(lh)>1
           [nodes, ~] = load_surface(lh{2});
           m.nodes=(m.nodes+nodes)/2; 
        end

        if ~isempty(rh) 
            idx_firstrh=size(m.nodes,1)+1;

            [nodes, triangles] = load_surface(rh{1});
            m.nodes = [m.nodes; nodes];
            m.triangles= [m.triangles; triangles+idx_firstrh-1];
            m.triangle_regions= [m.triangle_regions; 2*ones(size(triangles,1),1)];
        end
        if length(rh)>1
           [nodes, ~] = load_surface(rh{2});
           m.nodes(idx_firstrh:end,:)=(m.nodes(idx_firstrh:end,:)+nodes)/2; 
        end
        m.nodes=double(m.nodes);
        m.triangles=double(m.triangles);
    elseif (strcmp(model_creation_tool,'headreco') | strcmp(model_creation_tool,'mri2mesh'))
        m = mesh_load_fssurf(m2m_folder);
    end
    GM_mid_centers = (m.nodes(m.triangles(:,1),:) + m.nodes(m.triangles(:,2),:) + m.nodes(m.triangles(:,3),:))/3;
    GM_mid_centers = GM_mid_centers/1000;
end

function [nodes, triangles] = load_surface(fname)
    [~,~,extHlp] = fileparts(fname);
    if strcmpi(extHlp,'.gii') % load gifti
        s=gifti(fname);
        nodes=s.vertices;
        triangles=s.faces;
    else
        [nodes, triangles] = read_surf(fname);
        triangles=triangles+1; % FS indexing starts at 0
    end
end
