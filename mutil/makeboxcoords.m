function makeboxcoords(file_geom)
%==========================================================================
% When DIABOX adds air-sea fluxes to a box in your inverse model, it will 
% need to know if a gridded value from ERA fields is inside the box. To 
% allow for this, you must create box coordinate files, one for each box. 
% This routine allows easy creation of that polygon, by mouse-clicking on a 
% map of the sections and continents. 
%
% The program prompts the user for the box number, and saves the polygon 
% lons and lats in file "boxcoord/boxcoords*.mat" where * is the box number.
%
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
if (exist('file_geom')) 
   eval(file_geom);
else
   warning('Geometry file not found');
end

%% Plot map
clf;

figure('PaperType','A4','PaperUnits','centimeters',...
    'InvertHardcopy','off','Color',[1 1 1],'PaperPosition',[1 1 28.7 20],...
    'Position',[20 225 1300 705]);

Lon=[];Lat=[];
for i=1:size(geometry,2);
  eval(['load ',sectfiles(i,:),' lon lat']);
  % create a row vector.
  if size(lon,1) ~=1 % if lon and lat are column vectors
      Lon=[Lon,lon'];Lat=[Lat,lat'];
  else % if not
      Lon=[Lon,lon];Lat=[Lat,lat];
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


% Add box numbers in middle of each box.
for j=1:size(geometry,1);	% loop through boxes
  Lon=[];Lat=[];
  for i=1:size(geometry,2);	% accumulate lat, lon
    if (abs(geometry(j,i)));
      eval(['load ',sectfiles(i,:),' lon lat']);
      % create a row vector.
      if size(lon,1) ~=1 % if lon and lat are column vectors
          Lon=[Lon,lon'];Lat=[Lat,lat'];
      else % if not
          Lon=[Lon,lon];Lat=[Lat,lat];
      end

    end
  end
  Lon=mean(Lon);Lat=mean(Lat);
  text(Lon,Lat,num2str(j),'fontsize',15,'horizontalalignment','center', ...
	'verticalalignment','middle','color',colr(j,:));
end


%% Plot rivers
load dir_loc.mat riverdir;
eval(['load ',riverdir,'river_runoff.mat;']);
plot(lon_river,lat_river,'*b','markersize',6);

%% Create the polygons
zoom on;
s=input('Zoom to the desired box ("enter" when done)','s');
zoom off;
dim_lon = get(gca,'Xlim');
dim_lat = get(gca,'Ylim');

disp('Mouse-click to define box (Click outside the map when done when done)');
lat = []; lon = [];

count = 1;
while count == 1
    [xx,yy]=ginput(1);
    if (xx>dim_lon(2) | xx<dim_lon(1) | yy>dim_lat(2) | yy < dim_lat(1))
        count=0
    else
    lat = [lat yy]; 
    lon = [lon xx];
    count = 1;
    end
    clear xx yy s
end
plot([lon lon(1)],[lat lat(1)],'b');
hold off;

%% Save the data
boxnumber = input('Box number? ');
load dir_loc.mat boxcoord_dir
eval(['save ' boxcoord_dir 'boxcoords',num2str(boxnumber),' lon lat']);
return
