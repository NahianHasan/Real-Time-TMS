function [Eout,Hout,rcoil,kcoil,X,Y,Z,coil_thickness,min_coil,del_z] = volumetric_grid(p,del,rcoil,kcoil)
    del_z = max(p(3,:)) - min(rcoil(:,3));
        
    %Create a grid
    coil_thickness = abs(max(rcoil(:,3)) - min(rcoil(:,3)));
    min_coil = min(rcoil(:,3));
    del_scalp = [max(p(1,:))-min(p(1,:)),max(p(2,:))-min(p(2,:)),max(p(3,:))-min(p(3,:))];
    del_coil = [max(rcoil(:,1))-min(rcoil(:,1)),max(rcoil(:,2))-min(rcoil(:,2)),max(rcoil(:,3))-min(rcoil(:,3))];
    grid_extension = max(del_scalp);
    
    %[X,Y,Z] = ndgrid(min(p(1,:))-grid_extension:0.001:max(p(1,:))+grid_extension,min(p(2,:))-grid_extension:0.001:max(p(2,:))+grid_extension,min(p(3,:))-grid_extension:0.001:min(rcoil(:,3)) - (th_hair/2));
    [X,Y,Z] = ndgrid(-1.1*grid_extension:del:1.1*grid_extension,-1.1*grid_extension:del:1.1*grid_extension,-1.5*grid_extension:del:max(rcoil(:,3))+2*del);
    grid_loc = [reshape(X,[],1),reshape(Y,[],1),reshape(Z,[],1)];
    
    %Compute the field at grid points
    Ncoil=numel(rcoil)/3;
    Ns=numel(grid_loc)/3;
    disp('Calculating E Primary');
    [Eout]=computeEprimary(rcoil',kcoil',Ncoil,grid_loc',Ns);
    disp('Calculating H Primary');
    [Hout]=computeHprimary(rcoil',kcoil',Ncoil,grid_loc',Ns);
end