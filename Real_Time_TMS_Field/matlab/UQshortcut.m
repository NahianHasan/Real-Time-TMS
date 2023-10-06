function [Eaux_expected,Eaux_stddev]=UQshortcut(UQkit,posuq,anguq)

 Eaux_expected=UQkit.Eaux;
 Eaux_stddev=UQkit.Eaux;
 E_aux2=UQkit.Eaux.^2;
 for i=1:numel(UQkit.areap)
     measure=(UQkit.areap(:).*(UQkit.phi(:,i)<posuq));
     dmeasure=sum(measure);
     Eaux_expected(:,i)=UQkit.Eaux*measure;
     Eaux_expected(:,i)=Eaux_expected(:,i)/dmeasure;
     Eaux_stddev(:,i)=(E_aux2*measure)/dmeasure;
 end
Eaux_expected=Eaux_expected(:,UQkit.p);
Eaux_stddev=Eaux_stddev(:,UQkit.p);
if anguq~=0
a=ones([2*anguq+1 1])/(2*anguq);
a(1)=a(1)/2;a(end)=a(end)/2;
for i=1:numel(UQkit.p)
Eaux_expected([360-anguq+1:end 1:end-anguq],i) = cconv(a,Eaux_expected(:,i),360);
Eaux_stddev([360-anguq+1:end 1:end-anguq],i) = cconv(a,Eaux_stddev(:,i),360);
end
end
 Eaux_stddev=sqrt(Eaux_stddev-Eaux_expected.^2);
