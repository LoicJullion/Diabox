function plotsections(geometry,sectfiles)
%==========================================================================
% function plotsections(sm,cart)
%
% Plot all the sections in file "geometry file",
% adding labels giving the section names.
% AUTHOR: Loic Jullion (2015) based on Rick Lumpkin's file
%-----------------
% Licence:
% This file is licensed under the Creative Commons Attribution-Share 
% Alike 4.0 International license. 	
%
%     You are free:
% 
%         to share ? to copy, distribute and transmit the work
%         to remix ? to adapt the work
% 
%     Under the following conditions:
% 
%         attribution ? You must attribute the work in the manner specified 
%                       by the author or licensor (but not in any way that 
%                       suggests that they endorse you or your use of the 
%                       work).
%         share alike ? If you alter, transform, or build upon this work, 
%                       you may distribute the resulting work only under 
%                       the same or similar license to this one.
%-----------------
%
%==========================================================================
col='b';    % color

sm=1;
cart=1; 

Lon=[];Lat=[];
for i=1:size(geometry,2);
  eval(['load ',sectfiles(i,:),' lon lat']);
  % create a row vector.
  if size(lon,1) ~=1 % if lon and lat are column vectors
      Lon=[Lon,lon'];Lat=[Lat,lat'];
  else % if not
      Lon=[Lon,lon'];Lat=[Lat,lat'];
  end
  q=find(Lon<-180);Lon(q)=Lon(q)+360;
  q=find(Lon>180);Lon(q)=Lon(q)-360;
end;

hld=ishold;
if ~hld clf; end;

minlon = min([-100 min(Lon)-10]);

makemap([min(min(Lon),-100) max(Lon)],[min(Lat)-10 max(Lat)+10]);

hold on;
for i=1:size(geometry,2);	% loop through sections
  eval(['load ',sectfiles(i,:),' lon lat']);

  % detect sections wrapping around dateline
  q=find(abs(diff(lon))>180);
  while(length(q))
    lon=[lon(1:q(1)),NaN,lon(q(1)+1:length(lon))];
    lat=[lat(1:q(1)),NaN,lat(q(1)+1:length(lat))];
    q=find(abs(diff(lon))>180);
  end;
  if length(lon) 
     eval(['plot(lon,lat,''',col,''',''LineWidth'',2)']); 
  end;
  eval(['plot(lon(1),lat(1),''.'',''markersize'',10);']);
  % define median positions on section
  lon=lon+1e-9*rand(size(lon));
  lon2=lon(isfinite(lon));
  if (mod(length(lon2),2))
     lon2=lon2(1:length(lon2)-1);
  end;
  ind=find(lon2>=nanmedian(lon2));
  mlon=lon2(ind(1));mlat=lat(ind(1));

  name=sectfiles(i,:);
  q=find(name=='/');
  name(1:max(q))=[];
  q=find(name=='.');
  name(min(q):length(name))=[];
  text(mlon,mlat,name,'fontsize',14,'horizontalalignment','center');
end

if (~hld) hold off; end;

return;
