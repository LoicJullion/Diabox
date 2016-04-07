function makemap(lon,lat)
%==========================================================================
% function makemap(lon,lat);
% 
% Make a map spanning coordinates lon, lat (ex: lon=[-90 0]; lat=[-90 60];)
% for subsequent geophysical data. This function uses a low res version of
% Gebco (2.0 deg). This is intended for checking that the geometry of the
% box is correct, so the quality of the map is basic.
% 
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
load dir_loc.mat bathy_dir;
eval(['load ' bathy_dir 'gebco_2deg.mat']);

minlon = lon(1); maxlon = lon(2);
minlat = lat(1); maxlat = lat(2);

%imagesc(LON,LAT(end:-1:1),flipud(HH))
imagesc(lon_gebco,lat_gebco,HH_2deg)
hold on
contour(lon_gebco,lat_gebco,HH_2deg,[0 0],'Color',[1 1 1],'LineWidth',1)
axis xy
axis([minlon maxlon minlat maxlat]); caxis([-6000 0])
 
cmap =colormap(gray(64));
cmap = cmap(end:-1:1,:);

colormap(cmap);
 
