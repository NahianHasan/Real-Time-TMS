function [Edir]=coil2ROI(coilposid,ang,rs,ks,Anor,teid,te2p,p,conductivity,FEMord)
    rs2(1,:)=rs(1,:)*cos(ang/180*pi)+rs(2,:)*sin(ang/180*pi);
    rs2(2,:)=-rs(1,:)*sin(ang/180*pi)+rs(2,:)*cos(ang/180*pi);
    rs2(3,:)=rs(3,:);
    rs2(4,:)=1;
    ks2(1,:)=ks(1,:)*cos(ang/180*pi)+ks(2,:)*sin(ang/180*pi);
    ks2(2,:)=-ks(1,:)*sin(ang/180*pi)+ks(2,:)*cos(ang/180*pi);
    ks2(3,:)=ks(3,:);
    ks2=Anor(1:3,1:3,coilposid)*ks2;
    rs2=Anor(1:3,1:4,coilposid)*rs2;

    %run direct code and compute average along that
    rv=(p(:,te2p(1,:))+p(:,te2p(2,:))+p(:,te2p(3,:))+p(:,te2p(4,:)))/4;
    [Edirect]=runcodeksTensor(te2p,p,conductivity,rs2,ks2,rv,FEMord);%generates Efield at ro1

    Edir3=Edirect(:,teid);
    Edir_r=reshape(Edir3,1,[]);
    Edir=Edir_r.';
end