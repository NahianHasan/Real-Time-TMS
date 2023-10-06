function [Js,Ks,rho,rv,jv]=ROItoScalp(te2p,p,conductivity,teid,that,rs,nhat,FEMord)

%%%%find volume currents
[rv,jv]=runcoderecipmd(te2p,p,conductivity,teid-1,that,FEMord);

%%%find equivalent surface currents on rs,nhat
[Eout]=computeEprimary(rv,jv,numel(rv)/3,rs,numel(rs)/3);
[Hout]=computeHprimary(rv,jv,numel(rv)/3,rs,numel(rs)/3);
Hout=Hout/(4*pi*10^-7);
Ks=-cross(nhat,Eout,1);
Js=cross(nhat,Hout,1);
rho=sum(nhat.*Eout,1);

end