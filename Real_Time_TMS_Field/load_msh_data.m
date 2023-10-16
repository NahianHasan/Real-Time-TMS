function [p,te2p,conductivity,reg,M] = load_msh_data(msh_file,varargin)
    if strcmpi(msh_file(end-2:end),'mat')
        M = load(msh_file);
    elseif strcmpi(msh_file(end-2:end),'msh')
        M = feval(varargin{1},msh_file);
    else
        disp("Provide the msh file either in '.msh' format or '.mat' format. For '.msh' file, also give the function name");
        return;
    end
    p = double(M.nodes');
    p = p/1000;%conversion to m scale
    te2p = double(M.tetrahedra');
    reg = double(M.tetrahedron_regions);
    conductivity = reg;
    conductivity(conductivity==1) = 0.1260;
    conductivity(conductivity==2) = 0.2750;
    conductivity(conductivity==3) = 1.654;
    conductivity(conductivity==4) = 0.01;
    conductivity(conductivity==5) = 0.465;
    conductivity(conductivity==6) = 0.5;
    conductivity(conductivity==7) = 0.008;
    conductivity(conductivity==8) = 0.025;
    conductivity(conductivity==9) = 0.6;
    conductivity(conductivity==10) = 0.16;
    
    M.nodes = p';
    M.conductivity = conductivity;
    M.scalp_tri = surftri(p',te2p');%select scalp triangular facets
    M.original_msh_nodes = p';
end