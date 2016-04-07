function makemprop(geomfile,climato)
%==========================================================================
% function makemprop(geomfile);
% 
% Script which finds the mean temperature and salinity of a level surface 
% for an inverse model. It loads the relevant data from a given climatology 
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

% Load climatology. For now, only WGHC.
if strcmp(climato,'WGHC')==1
    % Load WGHC climatology
    % Available variables
    % LON LAT ZAX BOT_DEP PRES TEMP TPOTEN SALINITY OXYGEN SILICATE
    % NITRATE PHOSPH GAMMAN SIG0 SIG2 SIG4
    [lon_wghc,lat_wghc,ptemp,sal,gamma,press,bot_dpth] = ...
        read_wghc('LON','LAT','TPOTEN','SALINITY','GAMMAN','PRES',...
                  'BOT_DEP');
end

load dir_loc.mat boxcoord_dir diapvel_dir

% Create matrices containing the lat and lon;
latmat_wghc = repmat(lat_wghc,1,length(lon_wghc));
lonmat_wghc = repmat(lon_wghc',length(lat_wghc),1);

%% Map onto isopycnals
% define surfaces here
nsurfs = length(glevels);

%... find temperature and salinity values at each gamma
surf_ptemp = NaN*ones(length(lat_wghc),length(lon_wghc),length(glevels));
surf_salt  = surf_ptemp;
surf_press = surf_ptemp;

fprintf('Interpolating to find S, T, P on gamma surfaces.\n');
for ilat = 1:length(lat_wghc)
    for ilon = 1:length(lon_wghc);
        
        gp = squeeze(gamma(ilon,ilat,:));
        sp = squeeze(sal(ilon,ilat,:));
        tp = squeeze(sal(ilon,ilat,:));
        pp = squeeze(press(ilon,ilat,:));
        
        % NaN out any point where one of the variable is missing.
        q=find(isnan(gp+sp+tp));
        gp(q)=NaN;sp(q)=NaN;tp(q)=NaN;pp(q)=NaN;
        
        [surf_salt(ilat,ilon,:)]  = map2isopycnal(sp,gp,glevels);
        [surf_ptemp(ilat,ilon,:)] = map2isopycnal(tp,gp,glevels);
        [surf_press(ilat,ilon,:)] = map2isopycnal(pp,gp,glevels);
        
    end 
end
fprintf('Done \n');
fprintf('\n');

%% Loop over boxes
for boxnumber=1:nboxes;
    %--------------------------------
    % load lon, lat polygon defining box
    eval(['load ' boxcoord_dir 'boxcoords' num2str(boxnumber)]);

    % Calculate the area of the box 
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

    % calculate the area of the box in km^2
     total_area = polyarea(cumx,cumy);

    % convert to m^2
     total_area = total_area * 1e6;

    %-------------------------------- 
    % Find the grid points in the box
    % Look for the grid points of the WOCE climatology that are in the
    % polygon defined by the box.

    mark = inpolygon(lonmat_wghc,latmat_wghc,lon,lat);
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
    dx = 0.5*pi/180*re.*cos(latmat_wghc(box_point)*pi/180);
    d_area = sum(abs(dx.*dy));

    %--------------------------------
    % calculate areal averages of properties
    % Assign memory
    area    = NaN*ones(size(glevels)); % area of each layer
    msalt   = area;   mptemp = area;   mpres = area; % mean property
    stdsalt = area; stdptemp = area; stdpres = area; % stand. dev.
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
        plot(lonmat_wghc,latmat_wghc,'b.');
        hold on;
        plot(lon,lat,'k'); % Plot the box.
        % Plot the points within the box
        plot(lonmat_wghc(box_point),latmat_wghc(box_point),'r.')
        % Plot the points where the layer is present.
        plot(lonmat_wghc(box_point(layer_point)),latmat_wghc(box_point(layer_point)),'g.')
        subplot(3,1,2); % Now plot the corresponding map of temperature
        pcolor(lon_wghc,lat_wghc,squeeze(surf_ptemp(:,:,ilayer))); shading flat;
        hold on
        plot(lon,lat,'k');
        subplot(3,1,3); % And salinity
        pcolor(lon_wghc,lat_wghc,squeeze(surf_salt(:,:,ilayer))); shading flat;
        hold on
        plot(lon,lat,'k');

        % Make sure that the green points in subplot (3,1,1) correspond to the
        % points where temperature and salinity data exist.

        % Now calculate dx for the points within the box that contain the layer                                        
        dx     = 0.5*pi/180*re*cos(latmat_wghc(box_point(layer_point))*pi/180);
        area(ilayer) = sum(abs(dx.*dy)); % Calculate the area of the layer 
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

        % normalise so that maximum area within box does not exceed total area
        if max(area(ilayer))>total_area
           area(ilayer)=area(ilayer)*total_area/max(area);
        end

      if (ilayer==1)
        lwidth=surf_press(:,:,1);	% surface to level 1
      else
        lwidth=surf_press(:,:,ilayer)-surf_press(:,:,ilayer-1);
      end
      layer_width(ilayer)=nanmean(lwidth(find(lwidth>0)));
     %  
    end
    eval(['save ' diapvel_dir 'mprop_diap' num2str(boxnumber) '.mat mptemp msalt mpres stdsalt stdptemp stdpres layer_width area total_area ']);

end

return
%--------------------------------------------------------------------------
% OLD BITS OF CODE FROM RICK.

% lwidth=bot_dpth_box-reshape( ...
% 	surf_press(length(glevels),:,:) , ...
% 	[length(lat_sac) length(lon_sac)]);
% q=find(lwidth<=0);lwidth(q)=NaN;
% layer_width(length(glevels)+1)= ...
% 	nanmean(lwidth(find(lwidth>0)));

% % define surf_press2: equal to surf_press, but
% % NaNs replaced by 0 (surface too light to exist)
% % or bottom_sac value (surface is too dense).
% surf_press2=surf_press;
% for i=1:length(lat_sac);
%   for j=1:length(lon_sac);
%     q=find(isnan(surf_press2(:,i,j)));
%     if (length(q))
%       qq=q(find(q<length(glevels)/2));
%       if (length(qq)) surf_press2(qq,i,j)=0; end;
%       qq=q(find(q>length(glevels)/2));
%       if (length(qq)) surf_press2(qq,i,j)=bottom_sac(i,j);end;
%     end;
%   end;
% end;

% 
% % add finite (but small) areas for deep
% % and shallow layers, if they don't appear
% % in the SAC climatology.  For example,
% % dense layers exist in the Denmark Strait
% % section but not in the 1 deg. gridded
% % North Atlantic SAC fields.  A finite area
% % allows diapycnal fluxes in the inverse model.
% origarea=area;
% q=find(isnan(area) & defined_levels);
% if (length(q))
%   if (min(q) == min(find(defined_levels)))
%     qq=q(find(q<length(glevels)/2));
%     ind=min(qq):max(qq)+3;
%     x=ind;y=flipud(area(ind));
%     y2=fitdecay(x,y);
%     area(ind)=flipud(y2);
%     msalt(1:max(qq))=msalt(max(qq)+1);
%     mtemp(1:max(qq))=mtemp(max(qq)+1);
%     mpres(1:max(qq))=0;
%     stdsalt(1:max(qq))=stdsalt(max(qq)+1);
%     stdtemp(1:max(qq))=stdtemp(max(qq)+1);
%     stdpres(1:max(qq))=stdpres(max(qq)+1);
%   end;
%   if (max(q) == max(find(defined_levels)))
%     qq=q(find(q>length(glayers)/2));
%     ind=min(qq)-3:max(qq);
%     x=[ind,last(ind)+1];y=[area(ind);0];
%     y2=fitdecay(x,y);
%     area(ind)=y2(1:length(ind));
%     msalt(qq)=msalt(qq(1)-1);
%     mtemp(qq)=mtemp(qq(1)-1);
%     mpres(qq)=mpres(qq(1)-1);
%     stdsalt(qq)=stdsalt(qq(1)-1);
%     stdtemp(qq)=stdtemp(qq(1)-1);
%     stdpres(qq)=stdpres(qq(1)-1);
%   end;
%   if (min(q)>min(find(defined_levels)) & ...
% 	 max(q)<max(find(defined_levels)))
%     qq=find(isfinite(area));
%     area(q)=interp1(glevels(qq),area(qq), ...
% 	glevels(q),'linear');
%     msalt(q)=interp1(glevels(qq),msalt(qq), ...
% 	glevels(q),'linear');
%     mtemp(q)=interp1(glevels(qq),mtemp(qq), ...
% 	glevels(q),'linear');
%     mpres(q)=interp1(glevels(qq),mpres(qq), ...
% 	glevels(q),'linear');
%     stdsalt(q)=interp1(glevels(qq),stdsalt(qq), ...
% 	glevels(q),'linear');
%     stdtemp(q)=interp1(glevels(qq),stdtemp(qq), ...
%         glevels(q),'linear');
%     stdpres(q)=interp1(glevels(qq),stdpres(qq), ...
%         glevels(q),'linear');
%   end;
% end;
% 
% % now fix layer widths to be consistent with
% % the changes in level areas, assuming conservation
% % of volume and equipartition of layer width.
% origlayer_width=layer_width;
% w=interp1(glayers,layer_width,glevels,'linear');
% q=find(area~=origarea & isfinite(area));
% if (length(q))
%   qq=q(find(q<length(glevels)/2));
%   if (length(qq))
%     V=nansum(w(qq).*origarea(qq));
%     w(qq)=V./nansum(area(qq));
%   end;
%   qq=q(find(q>length(glevels)/2));
%   if (length(qq))
%     V=nansum(w(qq).*origarea(qq));
%     w(qq)=V./nansum(area(qq));
%   end;
% end;
% nlayer_width=[w(1);w(1:length(w)-1)+diff(w)/2;last(w)];
% q=find(isnan(layer_width));
% layer_width(q)=nlayer_width(q);
% 
% layer_area=[area(1); ...
% 	area(1:length(area)-1)+diff(area)/2; ...
% 	last(area)];
% 
% subplot(211);
% plot(origarea,'r'); hold on;
% plot([min(find(defined_levels)) ...
%   min(find(defined_levels))],[0 max(area)]);
% plot([max(find(defined_levels)) ...
%   max(find(defined_levels))],[0 max(area)]);
% plot(area,'.k');
% hold off;
% title('Level areas (corrected: points)')
% subplot(212);
% plot(origlayer_width,'r');hold on;
% plot(layer_width,'.k');hold off;
% title('Layer widths (corrected: points)');
% pause(5);
% %disp('Paused; enter "dbcont" to continue.');
% %keyboard;

%... save layer averages

%filename = input('box number (for filename suffix): ','s');
% filename=num2str(boxnumber);
% modelname='atl';
% 
% svar = [' msalt_',modelname,filename];
% tvar = [' mtemp_',modelname,filename]; 
% avar = [' area_',modelname,filename];
% pvar = [' mpres_',modelname,filename];
% svar2= [' stdsalt_',modelname,filename];
% tvar2= [' stdtemp_',modelname,filename];
% pvar2= [' stdpres_',modelname,filename];
% zvar2= [' layer_width_',modelname,filename];
% zvar3= [' layer_area_',modelname,filename];
% 
% eval([svar '=msalt;']);
% eval([tvar '=mtemp;']);
% eval([avar '=area;']);
% eval([pvar '=mpres;']);
% eval([svar2 '=stdsalt;']);
% eval([tvar2 '=stdtemp;']);
% eval([pvar2 '=stdpres;']);
% eval([zvar2 '=layer_width;']);
% eval([zvar3 '=layer_area;']);
% 
% eval(['save mprop_',modelname,filename svar tvar avar pvar ...
% 	svar2 tvar2 pvar2 zvar2 zvar3 ' lat_sac lon_sac ']);
% 
% % plot results
% subplot(234);
% plot(mtemp,-mpres,mtemp+stdtemp,-mpres,'r', ...
% 	mtemp-stdtemp,-mpres,'r');title('Temp')
% subplot(235);
% plot(msalt,-mpres,msalt+stdsalt,-mpres,'r', ...
% 	msalt-stdsalt,-mpres,'r');title('Salt');
% subplot(236);
% plot(area/1e12,-mpres);title('Area');
% 
% pause(1);
