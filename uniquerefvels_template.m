% uniquerefvels.m
%
% Define a set of inhomogeneous initial reference 
% velocities and a priori magnitudes of corrections for
% particular sections.

% directory containing dobox-prepared section data
% load defined levels, layers
glev;

load dir_loc.mat refvel_dir
%--------------------------------------------------------------------------
%% a08 (11 deg. S).
%
% Initial reference levels from Speer, Holfort,
%	Reynaud and Siedler (1996).
sect='a08.mat';
i=sectname2sectnum(sect,sectfiles);
if (length(i))

  refp=3900;
  eval(['load ',sect_dir,'/a08_raw bin_press surf_press']);
  d=(refp-nanmean(surf_press')').^2;
  reflevel(i)=find(d==min(d));

  Vmag=defVmag(sectfiles(i,:),reflevel(i),.01,100,.3);
  RV=rv(i)*ones(length(Vmag),1)';

  % shift reference level to isobaric surfaces
  % (ASSUMES BOTTOM REFERENCE LEVEL)
  reflevel(i)=length(glevels);  % must be set to bottom
                                % for isobaric code to work!
  eval(['load ',sect_dir,'/a08_raw bin_press lon binPair_vel']);              
  lon=lon(1:length(lon)-1)+diff(lon)/2;
  for j=1:length(lon);
    q=find(isfinite(binPair_vel(:,j)));
    binPair_vel(:,j)=binPair_vel(:,j)- ...
	binPair_vel(last(q),j);
  end;

  % shift ref level to 1100m west of 35.5W
  q=find(bin_press==1100);
  q2=find(lon<=-35.5);
  RV(q2)=RV(q2)-binPair_vel(q,q2);
  q=find(isnan(RV(q2)));RV(q2(q))=0;
  %RV(3)=-.2886;

  % shift ref level to 3800m between 35.5W and 17.67W
  % (Brazil Basin)  (stations 4:6 stay bottom-referenced)
  q=find(bin_press==3800);
  q2=find(lon>=-35.5 & lon<=-17.67);
  RV(q2)=RV(q2)-binPair_vel(q,q2);
  q=find(isnan(RV(q2)));RV(q2(q))=0;

  % shift ref level to 2400m between 17.67W and 10.17W (MAR)
  q=find(bin_press==2400);
  q2=find(lon>=-17.67 & lon<=-10.17);
  RV(q2)=RV(q2)-binPair_vel(q,q2);
  q=find(isnan(RV(q2)));RV(q2(q))=0;

  % shift ref level to 4000m east of 10.17W (Angola Basin)
  q=find(bin_press==4000);
  q2=find(lon>=-10.17);
  RV(q2)=RV(q2)-binPair_vel(q,q2);
  q=find(isnan(RV(q2)));RV(q2(q))=0;

  % uniform shift to balance Ekman (9.71Sv south)
  RV=RV-1.23e-4;

  clear binPair_vel lon bin_press q q2

  rvfilename=[refvel_dir,'/refvel_',sect];
  eval(['save ',rvfilename,' RV Vmag']);

end;

%--------------------------------------------------------------------------
%% a05fs (Florida Strait).
sect='a05fs.mat';
i=sectname2sectnum(sect,sectfiles);
if (length(i))

  reflevel(i)=length(glevels); %bottom

  Vmag=defVmag(sectfiles(i,:),reflevel(i),.3,-99,-99);
  eval(['load ',sectfiles(i,:),' binPair_vel']);
  Vmag=.15*ones(1,size(binPair_vel,2));

  RV=rv(i)*ones(length(Vmag),1)';
  RV(1:length(RV))=.285;  % 31 Sv through Florida Strait
  for bl=1:length(Vmag);
    Vmag(bl)=.15;
  end;

  rvfilename=[refvel_dir,'/refvel_',sect];
  eval(['save ',rvfilename,' RV Vmag']);
end;

%--------------------------------------------------------------------------
%% a05int (24.5 deg. N).
%   Increase the magnitude of reference
%   velocities against the western slope
sect='a05int.mat';
i=sectname2sectnum(sect,sectfiles);
if (length(i))

  % 3000 dbar ref level (Rintoul and Wunsch 1991)
  refp=3000;
  eval(['load ',sect_dir,'/a05int_raw bin_press surf_press']);
  d=(refp-nanmean(surf_press')').^2;
  reflevel(i)=find(d==min(d));
  Vmag=defVmag(sectfiles(i,:),reflevel(i),.01,100,.3);
  RV=rv(i)*ones(length(Vmag),1)';

  % shift RV to zero the thermal wind profile at 3000dbar
  reflevel(i)=length(glevels); 	% must be set to bottom
				% for isobaric code to work!
  eval(['load ',sect_dir,'/a05int_raw binPair_vel lon bin_press bottom']);
  lon=lon(1:length(lon)-1)+diff(lon)/2;
  bottom=[bottom(1:length(bottom)-1),bottom(2:length(bottom))];
  bottom=nanmin(bottom')';
  for j=1:length(lon);
    q=find(isfinite(binPair_vel(:,j)));
    binPair_vel(:,j)=binPair_vel(:,j)- ...
        binPair_vel(last(q),j);
  end;
  q=find(bin_press==refp);
  q2=find(bottom>refp);
  RV(q2)=RV(q2)-binPair_vel(q,q2);
  q=find(isnan(RV(q2)));RV(q2(q))=0;

  % balance Ekman and FS flow: want 35.3 Sv south:
  %	4.3Sv N Ekman + 31.0 Sv N in Florida Strait
  RV=RV-8.7e-4;

  rvfilename=[refvel_dir,'/refvel_',sect];
  eval(['save ',rvfilename,' RV Vmag']);
end;

%--------------------------------------------------------------------------
%% ar16d : Gulf of Cadiz
sect='ar16d.mat';
i=sectname2sectnum(sect,sectfiles);
if (length(i))
  reflevel(i)=length(glevels);	% bottom
  Vmag=defVmag(sectfiles(i,:),reflevel(i),.1,-99,-99);
  RV=rv(i)*ones(length(Vmag),1)';

  eval(['load ',sect_dir,'/ar16d_raw lat']);
  lat=lat(1:length(lat)-1)+diff(lat)/2;

  % shift to cancel net flow at sub-NACW densities
  % (<7 degC; gamma>=27.94)
  q=find(lat>33.85 & lat<36.2); RV(q)=RV(q)-5e-3;
  q=find(lat>34 & lat<34.5); RV(q)=5e-3;
  % 2 Sv outflow in layers 15-20, met by magnifying
  % westward flow against Spain
  q=find(lat>36.2 & lat<36.75); RV(q)=.02;
  % eastward flow of surface water against Africa
  % to balance net Ekman (0.15Sv west)
  q=find(lat<33.85); RV(q)=-1.6e-1;

  rvfilename=[refvel_dir,'/refvel_',sect];
  eval(['save ',rvfilename,' RV Vmag']);
end;

%--------------------------------------------------------------------------
%% ar19a: 48N

includeNAC=1;   % set to 0 for LNM at neutral density 27.879

sect='ar19a.mat';
i=sectname2sectnum(sect,sectfiles);
if (length(i))
  % adjustments to reference velocities from
  % Meinen, 2001: DSR I 48, 1553-1580.
  q=find(glevels>=27.879); reflevel(i)=q(1);
  if ~includeNAC
    Vmag=defVmag(sectfiles(i,:),reflevel(i),.02,100,.3);
    RV=zeros(length(Vmag),1)';
  else

  eval(['load ',sect_dir, ...
          '/ar19a lon bin_press bin_temp lon lat bottom;'])
  defVmag_48n;
  RV=zeros(length(Vmag),1)';

  
  mann=-43.85;  % longitude of Mann Eddy center
  bottom=bottom(1:length(lon)-1)+diff(bottom)/2;
  
  % Identify the core of the NAC (where the 10dgC
  % isotherm crosses 450dbar)
  temp=bin_temp(find(bin_press==450),:);
  lonhr=min(lon):.01:max(lon);
  lathr=interp1(lon,lat,lonhr,'linear');
  temp=interp1(lon,temp,lonhr,'linear');
  lon_0=lonhr(min(find(temp>10)));
  lat_0=lathr(min(find(temp>10)));
  % convert lon, lat to pair coordinates
  lon=lon(1:length(lon)-1)+diff(lon)/2;
  lat=lat(1:length(lat)-1)+diff(lat)/2;
  
  % calculate distance (negative=west)
  % between station pairs and NAC core.
  dist=NaN*ones(size(lon));
  for jj=1:length(dist);
    dist(jj)=sw_dist([lat(jj) lat_0], ...
        [lon(jj) lon_0],'km');
    if (lon(jj)<lon_0) dist(jj)=-dist(jj); end;
  end;
  
  % load mean reference velocity on 27.879 surface
  %  (calculated by "calcstreammean.m" in NACpies directory)
  load ~/Work/matlab/dobox/new_dobox_mfiles/48n_refvel d_s v27879
  % rescale so that max(d_s) is at mann eddy center 
  % for this repeat (instead of +200km)
  q=find(d_s>0);
  d_s(q)=d_s(q)/200*interp1(lon,dist,mann);
  % interpolate to dist for this repeat
  v2=interp1(d_s,v27879,dist,'linear');
  % extrapolate westward to shelf break (bottom==500)
  q=find(dist<0 & bottom'>500 & isnan(v2));
  v2(q)=v27879(1);
  % zero bottom ref level inshore of this
  q=find(dist<0 & isnan(v2)); v2(q)=0;
  % zero ref level east of mann eddy center
  v2(find(dist>interp1(lon,dist,mann)))=0;
  
  % match the time-mean net transport from -200km to Mann 
  % center (pairs 4:14):
  %   145 Sv northward - 30 Sv southward = 115 Sv
  %q=find(dist>-200 & lon<mann);
  v2=0.802*v2;
  RV=RV+v2;
  
  % East half of Mann eddy:  55Sv south.
  q=find(lon>=-43.85 & lon<=-41);  % pairs 15:18
  RV(q)=RV(q)-4.21e-2;
  
  % these adjustments give a net transport of 64 Sv north.
  % adjust in interior (west of eastern edge Mann eddy)
  % to get +2.46 Sv (-1*Ekman transport)
  q=19:length(lon); 
  RV(q)=RV(q)-6.59e-3;
  
  rvfilename=[refvel_dir,'/refvel_',sect];
  eval(['save ',rvfilename,' RV Vmag']);
  
 end; % if includeNAC
  
end;


clear sect i;
