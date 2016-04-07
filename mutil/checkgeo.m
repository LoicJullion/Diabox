function checkgeo(file_geom)
%==========================================================================
% checkgeo.m
%
% graphically check the geometry file for a box
% to make sure the sections and corresponding
% "geometry" values are correct.
% AUTHORS: Loic Jullion based on Rick Lumpkin's code.
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

figure('PaperType','A4','PaperUnits','centimeters',...
    'InvertHardcopy','off','Color',[1 1 1],'PaperPosition',[1 1 28.7 20],...
    'Position',[20 225 1300 705]);

clf;

if (exist('file_geom')) 
   eval(file_geom)

else
   warning('Geometry file not found');
end;

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

makemap([min(Lon)-10 max(Lon)+10],[min(Lat)-10 max(Lat)+10]);

colr=[[1 0 0];[0 0 1];[.14 .55 .14]; ...
	[.8 .8 0];[.8 0 .8];[.5 .5 .5]; ...
	[.87 .72 .53];[1 .5 .3]; ...
	[1 .4 .7];[1 .84 0]; ...
	[.94 .9 .55];[.13 .7 .67]; ...
	[0 0 .5];[.55 .27 .07]; ...
	[0 .5 .5];[1 .39 .28]];
[n,i]=sort(rand(size(colr,1),1));
colr=colr(i,:);
if (size(colr,1)<size(geometry,1))
  colr=[colr;rand( ...
	size(geometry,1)-size(colr,1),3)];
end;
colr=colr(1:size(geometry,1),:);

hold on;
plotsections(geometry,sectfiles);

for i=1:size(geometry,2);	% loop through sections
  eval(['load ',sectfiles(i,:),' lon lat']);
  theta=360/2/pi*atan2( median(diff(lat)),median(diff(lon)));

  % define median positions on section
  ind=find(lat==min(lat(find(lat>=nanmedian(lat)))));
  mlon=lon(ind(1)); mlat=lat(ind(1));

  for j=1:size(geometry,1);	% loop through boxes
    if (geometry(j,i)~=0)
      mf=.05*(range(Lon)+range(Lat))/2;
      theta2=theta*2*pi/360+geometry(j,i)*pi/2;
      theta2=theta2+.1*randn(1,1);
      P=quiver([mlon 0],[mlat 0], ...
	    mf*[cos(theta2) 0],mf*[sin(theta2) 0],0);
      set(P,'color',colr(j,:),'linewidth',1);
    end
  end
  plot(mlon,mlat,'.w','markersize',15);
  text(mlon,mlat,num2str(i),'fontsize',8, ...
	'horizontalalignment','center');
end;
% Add box numbers in middle of each box.
for j=1:size(geometry,1);	% loop through boxes
  Lon=[];Lat=[];
  for i=1:size(geometry,2);	% accumulate lat, lon
    if (abs(geometry(j,i)));
      eval(['load ',sectfiles(i,:),' lon lat']);
      Lon=[Lon,lon];Lat=[Lat,lat];
    end
  end
  Lon=mean(Lon);Lat=mean(Lat);
  text(Lon,Lat,num2str(j),'fontsize',15,'horizontalalignment','center', ...
	   'verticalalignment','middle', 'color',colr(j,:));
end
hold off;
disp('Arrows should point into boxes.');

titl=pwd;
title(['Geometry of ',titl],'fontsize',15);

clear Lat Lon P colr d i j lat lon mlat mlon
clear theta theta2 titl

zoom on;
