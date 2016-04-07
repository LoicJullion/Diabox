function g = gamman(s,t,p,lon,lat);

% function g = gamman(s,t,p,lon,lat);
%
% s=practical salinity.  
% t=IN SITU temperature (converted to potential temp below).
% p=pressure (dbar)
%
% Call gamma_n, after adjusting for lats and
% lons outside of the range of the gamma_n routine.
% Also, if "p" is a vector, convert it automatically
% to a matrix corresponding to s, t (length(p) by 
% length(lon)).
%
% Version 2: 16 July 2012: use gamma_GP_from_SP_pt.m
%
% test value: gamman(34.9,3.6,2000,-49,43)=27.9280
%   with version 2: 27.926

q=find(lon>180);lon(q)=lon(q)-360;
% Move points inside inland seas to be just outside
% those seas, to avoid places where gamma_n fails.
SLon=[-6 -6 14 43 34]; SLat=[30 41 49 37 30];
q=find(isinpoly(lon,lat,SLon,SLat));
if (length(q))
  disp('Moving data within Med.Sea to Gulf of Cadiz')
  lon(q)=-10;lat(q)=35;
end;
SLon=[9 9 24 35 23];SLat=[52 61 68 58 52];
q=find(isinpoly(lon,lat,SLon,SLat));
if (length(q))
  disp('Moving data within Baltic Sea to North Sea')
  lon(q)=7;lat(q)=57;
end;
SLon=[39 27 38 48];SLat=[30 30 9 16];
q=find(isinpoly(lon,lat,SLon,SLat));
if (length(q))
  disp('Moving data within Red Sea to Aden Gulf')
  lon(q)=46;lat(q)=12;
end;
SLon=[50 46 54 58];SLat=[33 29 21 30];
q=find(isinpoly(lon,lat,SLon,SLat));
if (length(q))
  disp('Moving data within Persian Gulf to Omen Gulf')
  lon(q)=58;lat(q)=25;
end;
SLon=[-71 -78 -98 -85];SLat=[69 47 59 73];
q=find(isinpoly(lon,lat,SLon,SLat));
if (length(q))
  disp('Moving data within Hudson Bay to Lawrence')
  lon(q)=-71;lat(q)=62;
end;
SLon=[115.1 115.1 115.5 115.5];SLat=[-8 -10 -10 -8];
q=find(isinpoly(lon,lat,SLon,SLat));
if (length(q))
  disp('Moving data against Indonesia to south')
  lon(q)=115.5;lat(q)=-10;
end;

% shift points in Arctic Ocean
q=find(lat>68 & lon>=16 & lon<=78); 
  lon(q)=15.5;lat(q)=70;
q=find(lat>68 & lon>=78 & lon<=140); lon(q)=140.5;
q=find(lat>65 & lon>=-136 & lon<=-106); 
  lon(q)=-136.5;lat(q)=80;
q=find(lat>72 & lon>=-106 & lon<=-76); 
  lon(q)=-75.5;lat(q)=86;
q=find(lat>70 & lon>=-160 & lon<=-158); lon(q)=-160.5;
q=find(lat>70 & lon>=-158 & lon<=-156); lon(q)=-155.5;

alat=lat;
q=find(lat<-79.999);
if (length(q));
  fprintf('Warning: lats below 80S considered 80S by gamma_n.\n');
  fprintf('  Minimum latitude: %f\n',min(lat(q)));
  alat(q)=-79.999;
end;
q=find(lat>63.999);
if (length(q));
  fprintf('Warning: lats above 64N considered 64N by gamma_n.\n');
  fprintf('  Maximum latitude: %f\n',max(lat(q)));
  alat(q)=63.999;
end;

alon=lon;
q=find(alon<0);
alon(q)=alon(q)+360;

if size(p,2)==1
  p=p*ones(1,length(lon));
end;

g = gamma_GP_from_SP_pt(s,sw_ptmp(s,t,p,0),p,alon,alat);

% regress onto potential density
% to fill the bad values

maxg=max(max(g));
q=find(g>0);
ming=min(min(g(q)));
q=find(g<0);
if (length(q));
  disp('gamma_n failed in some bins: filling by regression.');
end;
for i=1:length(q);
  fprintf('%d of %d\n',i,length(q));
  rho=sw_pden(s,t,p,p(q(i)));
  drho=.1;
  qq=find(rho>=rho(q(i))-drho & ...
	rho<=rho(q(i))+drho & g>0);
  if (length(find(g>0))<3) 
    disp('No good values of gamma for regression!')
    break;
  end;
  while (length(qq)<3)
    drho=drho+.1;
    qq=find(rho>=rho(q(i))-drho & ...
        rho<=rho(q(i))+drho & g>0);
  end;
  pf=polyfit(rho(qq),g(qq),1);
  g(q(i))=polyval(pf,rho(q(i)));
  if (g(q(i))>maxg) 
    disp('  interpolation out-of-range (maxg)')
   g(q(i))=maxg; 
  elseif (g(q(i))<ming)
      disp(' interpolation out-of-range (ming)')
      g(q(i))=ming;
  end;
end;

q=find(g<0);
if (length(q))
  disp('Warning: interpolation routine failed.');
  g2=g;
  g2(q)=NaN;
  g2=fillnans(g2,0);
  g2=g2.*g./g;
  g=g2;
end;

