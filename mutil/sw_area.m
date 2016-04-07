function [Area,dx,dy]=sw_area(lon,lat)
%==========================================================================
% function [Area,dx,dy]=sw_area(lon,lat);
% 
% Given vectors lon, lat (degrees), calculate matrix Area (units: m^2) of 
% dimension length(lon) by length(lat). 
%
% AUTHORS: Rick Lumpkin
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
if (length(lon)==size(lon,2))
  lon=lon';
end;
if (length(lat)==size(lat,2))
  lat=lat';
end;

dlon=diff(lon); dlat=diff(lat);

% define lon2, lat2 bracketing values of
% lon and lat.  length(lon2)=length(lon)+1
lon2=lon(1:length(lon)-1)+dlon/2;
lon2=[lon(1)-dlon(1)/2;lon2;last(lon)+last(dlon)/2];
lat2=lat(1:length(lat)-1)+dlat/2;
lat2=[lat(1)-dlat(1)/2;lat2;last(lat)+last(dlat)/2];

if (max(diff(lon)) < min(diff(lon))/.99999);
	% evenly-spaced longitudinal grid; fast calculation
    % .99999 added because 1/3 grid + truncation noise
    % caused max(diff(lon))/min(diff(lon))=1 + 1.7e-13

  Area=NaN*ones(length(lat),1);
  dx=Area;dy=Area;
  for i=1:length(lat);
    dx(i)=sw_dist([lat(i) lat(i)], ...
      [lon2(1) lon2(2)],'km')*1e3;
    dy(i)=sw_dist([lat2(i) lat2(i+1)], ...
      [lon(1) lon(1)],'km')*1e3;
  end;
  Area=dx.*dy;
  Area=Area*ones(1,length(lon));
  dx=dx*ones(1,length(lon));
  dy=dy*ones(1,length(lon));

else;	% uneven lon grid; loop through both
    
  Area=NaN*ones(length(lat),length(lon));
  dx=Area;dy=Area;
  for i=1:length(lat);
    for j=1:length(lon);
      dx(i,j)=sw_dist([lat(i) lat(i)], ...
	[lon2(j) lon2(j+1)],'km')*1e3;
      dy(i,j)=sw_dist([lat2(i) lat2(i+1)], ...
	[lon(j) lon(j)],'km')*1e3;
    end;
  end;
  Area=dx.*dy;

end

