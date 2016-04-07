function boxnum=getboxnum(lonp,latp,nboxes)
%==========================================================================
% function boxnum=getboxnum(lon,lat,nboxes);
%
%  given a lon and lat, find the box which contains that point and output
%  the box number.
%
% AUTHORS: Rick Lumpkin.
%
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
if ~exist('nboxes')
    geo
    nboxes=size(geometry,1);
end;

if lonp>180 lonp=lonp-360; end;
load dir_loc.mat boxcoord_dir
boxnum=[];
for bn=1:nboxes;
  eval(['load ' boxcoord_dir 'boxcoords',num2str(bn)]);
  if (isinpoly(lonp,latp,lon,lat))
    boxnum=bn;
    return;
  end;
end;
return