function makemprop(geomfile,climato)
%==========================================================================
% function makemprop(geomfile,climato);
% 
% Script which finds the mean temperature and salinity of a given  
% isopycnal surface. It loads the relevant data from a given climatology 
% (WHICH ONE?).
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
%% Load data
% load geometry
eval(geomfile);
nboxes=size(geometry,1);
% load isopycnals
glev;
% Radius of the earth (In m).
re=6370e3; 

load dir_loc.mat boxcoord_dir diapvel_dir

%% Load climatology. For now, only WGHC.
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% This part of the code needs to be updated
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if strcmp(climato,'WGHC')==1
    % Load WGHC climatology
    % Available variables
    % LON LAT ZAX BOT_DEP PRES TEMP TPOTEN SALINITY OXYGEN SILICATE
    % NITRATE PHOSPH GAMMAN SIG0 SIG2 SIG4
    [lon_clim,lat_clim,temp,ptemp,sal,gamma,press,bot_dpth] = ...
        read_wghc('LON','LAT','TEMP','TPOTEN','SALINITY','GAMMAN','PRES',...
                  'BOT_DEP');
              
    if isempty(gtype) == 1 
       % Neutral density. Already calculated in read_wghc
       
    else
       % To be modified depending on which density is being used 
       temp1 = reshape(temp,length(lat_clim)*length(lon_clim),size(temp,3));
       sal1 = reshape(sal,length(lat_clim)*length(lon_clim),size(temp,3));
       press1 = reshape(press,length(lat_clim)*length(lon_clim),size(temp,3));
       
       sig0    = temp1.*NaN;
       sig1000 = temp1.*NaN;
       
       sig0    = sw_pden(sal1,temp1,press1,0)-1000.;
       sig1000 = sw_pden(sal1,temp1,press1,1000)-1000.;
       
       sig0 = reshape(sig0,[length(lon_clim),length(lat_clim),size(temp,3)]);
       sig1000 = reshape(sig1000,[length(lon_clim),length(lat_clim),size(temp,3)]);
       
    end
end

% Create matrices containing the lat and lon;
latmat_clim = repmat(lat_clim,1,length(lon_clim));
lonmat_clim = repmat(lon_clim',length(lat_clim),1); 

%% Map onto isopycnals
% define surfaces here
nsurfs = length(glevels);

%... find temperature and salinity values at each gamma
surf_ptemp = NaN*ones(length(lat_clim),length(lon_clim),length(glevels));
surf_salt  = surf_ptemp;
surf_press = surf_ptemp;

fprintf('Interpolating to find S, T, P on gamma surfaces.\n');
for ilat = 1:length(lat_clim)
    for ilon = 1:length(lon_clim);
        if isempty(gtype) == 1
            gp = squeeze(gamma(ilon,ilat,:));
            sp = squeeze(sal(ilon,ilat,:));
            tp = squeeze(ptemp(ilon,ilat,:));
            pp = squeeze(press(ilon,ilat,:));

            % NaN out any point where one of the variable is missing.
            q=find(isnan(gp+sp+tp));
            gp(q)=NaN;sp(q)=NaN;tp(q)=NaN;pp(q)=NaN;

            [surf_salt(ilat,ilon,:)]  = map2isopycnal(sp,gp,glevels);
            [surf_ptemp(ilat,ilon,:)] = map2isopycnal(tp,gp,glevels);
            [surf_press(ilat,ilon,:)] = map2isopycnal(pp,gp,glevels);

            clear gp sp tp pp q
        else
            % This needs to be changed to reflect user's choices 
            s0    = squeeze(sig0(ilon,ilat,:));
            s1000 = squeeze(sig1000(ilon,ilat,:));
            sp    = squeeze(sal(ilon,ilat,:));
            tp    = squeeze(ptemp(ilon,ilat,:));
            pp    = squeeze(press(ilon,ilat,:));

            % NaN out any point where one of the variable is missing.
            q=find(isnan(s0+s1000+sp+tp));
            s0(q)=NaN;s1000(q)=NaN;sp(q)=NaN;tp(q)=NaN;pp(q)=NaN;

            [surf_salt_0(ilat,ilon,:)]  = map2isopycnal(sp,s0,glevels(1:2));
            [surf_ptemp_0(ilat,ilon,:)] = map2isopycnal(tp,s0,glevels(1:2));
            [surf_press_0(ilat,ilon,:)] = map2isopycnal(pp,s0,glevels(1:2));
            
            [surf_salt_1000(ilat,ilon,:)]  = map2isopycnal(sp,s1000,glevels(3));
            [surf_ptemp_1000(ilat,ilon,:)] = map2isopycnal(tp,s1000,glevels(3));
            [surf_press_1000(ilat,ilon,:)] = map2isopycnal(pp,s1000,glevels(3));
            
            % Concatenate the mapped properties together
            surf_salt  = cat(3,surf_salt_0,surf_salt_1000);
            surf_ptemp = cat(3,surf_ptemp_0,surf_ptemp_1000);
            surf_press = cat(3,surf_press_0,surf_press_1000);          
        end
    end 
end
fprintf('Done \n');
fprintf('\n');
%% Loop over boxes
for boxnumber=1:nboxes;
    %--------------------------------
    % load lon, lat polygon defining box
    eval(['load ' boxcoord_dir 'boxcoords' num2str(boxnumber)]);

    % Calculate the area_layer of the box 
    % calculate distance and angle between consecutive stations in the box
     [dist,phaseangle] = sw_dist(lat,lon,'km');

    % convert phaseangle to radians
     phaseangle = phaseangle *(2*pi/360);

    % convert these to increments of x and y
     x = dist.*(cos(phaseangle));
     y = dist.*(sin(phaseangle));

    % make the first point the origin
     cumx = [0 cumsum(x)];
     cumy = [0 cumsum(y)];

    % calculate the area_layer of the box in km^2
     total_area = polyarea(cumx,cumy);

    % convert to m^2
     total_area = total_area * 1e6;

    %-------------------------------- 
    % Find the grid points in the box
    % Look for the grid points of the WOCE climatology that are in the
    % polygon defined by the box.

    mark = inpolygon(lonmat_clim,latmat_clim,lon,lat);
    % Extract the index of the points within the box
    box_point = find(mark == 1 | mark==0.5);

    % The arc distance between two points sitting on the same meridian is:
    % dy = R * d Phi, where R is the Earth Radius and d Phi is the difference
    % in latitude.
    % The arc distance between two points on the same latitude is:
    % dx = R * cos(Phi) * d Lambda, where R is the Earth Radius, d Lambda is
    % the longitude difference and  cos(Phi) is the cosine of the latitude.
    %
    % Here we have a 0.5 degrees resolution:
    dy = 0.5*pi/180*re; 
    dx = 0.5*pi/180*re.*cos(latmat_clim(box_point)*pi/180);
    d_area = sum(abs(dx.*dy));

    %--------------------------------
    % calculate areal averages of properties
    % Assign memory
    area_layer    = NaN*ones(size(glevels)); % area_layer of each layer
    msalt   = area_layer;   mptemp = area_layer;   mpres = area_layer; % mean property
    stdsalt = area_layer; stdptemp = area_layer; stdpres = area_layer; % stand. dev.
    layer_width = NaN*ones(length(glevels)+1,1);
    disp('Calculating layer mean values.');
    for ilayer=1:length(glevels);

        ptempmap = squeeze(surf_ptemp(:,:,ilayer));
        salmap   = squeeze(surf_salt(:,:,ilayer));
        pressmap = squeeze(surf_press(:,:,ilayer));

        % Keep only the points in the box
        ptempmap  = ptempmap(box_point); 
        salmap    = salmap(box_point);
        pressmap  = pressmap(box_point);

        layer_point = find(isnan(ptempmap)==0); % Exclude the grid points where 
                                                % the layer interface is not
                                                % present                            
        % Plot some diagnostics
        figure;
        subplot(3,1,1);
        plot(lonmat_clim,latmat_clim,'b.');
        hold on;
        plot(lon,lat,'k'); % Plot the box.
        % Plot the points within the box
        plot(lonmat_clim(box_point),latmat_clim(box_point),'r.')
        % Plot the points where the layer is present.
        plot(lonmat_clim(box_point(layer_point)),latmat_clim(box_point(layer_point)),'g.')
        subplot(3,1,2); % Now plot the corresponding map of temperature
        pcolor(lon_clim,lat_clim,squeeze(surf_ptemp(:,:,ilayer))); shading flat;
        hold on
        plot(lon,lat,'k');
        subplot(3,1,3); % And salinity
        pcolor(lon_clim,lat_clim,squeeze(surf_salt(:,:,ilayer))); shading flat;
        hold on
        plot(lon,lat,'k');

        % Make sure that the green points in subplot (3,1,1) correspond to the
        % points where temperature and salinity data exist.

        % Now calculate dx for the points within the box that contain the layer                                        
        dx     = 0.5*pi/180*re*cos(latmat_clim(box_point(layer_point))*pi/180);
        area_layer(ilayer) = sum(abs(dx.*dy)); % Calculate the area_layer of the layer 
                                           % within the box (m^2).

        % Now calculate the averaged temperature and salinity for the layers
        % within the box.
        % Mean properties for each isopycnal
        mptemp(ilayer) = nanmean(ptempmap(layer_point));
        msalt(ilayer)  = nanmean(salmap(layer_point));
        mpres(ilayer)  = nanmean(pressmap(layer_point));
        % STD.
        stdptemp(ilayer) = nanstd(ptempmap(layer_point));
        stdsalt(ilayer)  = nanstd(salmap(layer_point));
        stdpres(ilayer)  = nanstd(pressmap(layer_point));

        % normalise so that maximum area_layer within box does not exceed total area_layer
        if max(area_layer(ilayer))>total_area
           area_layer(ilayer)=area_layer(ilayer)*total_area/max(area_layer);
        end

      if (ilayer==1)
        lwidth=surf_press(:,:,1);	% surface to level 1
      else
        lwidth=surf_press(:,:,ilayer)-surf_press(:,:,ilayer-1);
      end
      layer_width(ilayer)=nanmean(lwidth(find(lwidth>0)));
     %  
    end
    eval(['save ' diapvel_dir 'mprop_diap' num2str(boxnumber) '.mat mptemp msalt mpres stdsalt stdptemp stdpres layer_width area_layer total_area ']);

end

return
%--------------------------------------------------------------------------
