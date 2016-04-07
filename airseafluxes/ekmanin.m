function [EkmanInM,EkmanInH,EkmanInS]=ekmanin(lon_coads,lat_coads, ...
	Ue,Ve,Gamma,Temp,Salam,glevels,file_geom,boxnumber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [EkmanInM,EkmanInH,EkmanInS]=ekmanin(lon_coads,lat_coads, ...
%        Ue,Ve,Gamma,Temp,Salam,glevels);
%
% Given the 2D fields Ue and Ve (m^2/s) and the corresponding
% fields Gamma, Temp, Salam, calculate the net Ekman transports
% into the box as a function of gamma (in bins defined by glevels).
%
% Ekman velocities are rotated to be in an along- and across-section
% direction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% call box geometry file
eval(file_geom);
% isolate line of geometry matrix for this box
geometry=geometry(boxnumber,:);

% allocate memory for net Ekman transports
EkmanInM=zeros(length(glevels),size(sectfiles,1));
EkmanInH=EkmanInM; EkmanInS=EkmanInM;

% identify sections which border the box
bsect=find(geometry~=0);
for sectind=1:length(bsect);
  sect=bsect(sectind);

  % load section's station locations
  eval(['load ',sectfiles(sect,:),' lon lat']);
  % calculate position, separation distance 
  % and angle of each station pair
  lon2=lon(1:length(lon)-1) + diff(lon)/2;
  lat2=lat(1:length(lat)-1) + diff(lat)/2;
  for stn=1:length(lon2)
    if (lon2(stn)>max(lon_coads)) lon2(stn)=max(lon_coads);end;
    if (lon2(stn)<min(lon_coads)) lon2(stn)=min(lon_coads);end;
    if (lat2(stn)>max(lat_coads)) lat2(stn)=max(lat_coads);end;
    if (lat2(stn)<min(lat_coads)) lat2(stn)=min(lat_coads);end;
  end;
  [r,theta]=sw_dist(lat,lon,'km');
  r=r*1e3; theta=theta*2*pi/360;

  % interpolate u, v, and 10m values onto stn pair locations
  u=interp2(lon_coads,lat_coads,fillnans(Ue,0), ...
	lon2,lat2,'linear');
  v=interp2(lon_coads,lat_coads,fillnans(Ve,0), ...
	lon2,lat2,'linear');
  gamma_ek=interp2(lon_coads,lat_coads,fillnans(Gamma,0), ...
        lon2,lat2,'linear');
  temp_ek=interp2(lon_coads,lat_coads,fillnans(Temp,0), ...
        lon2,lat2,'linear');
  salam_ek=interp2(lon_coads,lat_coads,fillnans(Salam,0), ...
        lon2,lat2,'linear');

  % calculate EkIn: pair-by-pair Ekman transport
  % into box (at this point, units are m^2/s)
  EkIn=(-u.*sin(theta)+v.*cos(theta));

  % plot Ekman transports (2D field and interpolated values)
  quiver(lon_coads,lat_coads,Ue,Ve,0)
  hold on;
  plot(lon,lat,'.k')
  u2=EkIn.*cos(theta+pi/2);v2=EkIn.*sin(theta+pi/2);
  P=quiver(lon2,lat2,u2,v2,0);
  set(P,'color','r');
  set(gca,'dataaspectratio',[1 1 1]);
  hold off;
  axis([min(lon)-3 max(lon)+3 min(lat)-10 max(lat)+10]);
  title(['Ekman flux across ',sectfiles(sect,:)]);
  pause(.1);

  % convert EkIn to m^3/s, including value
  % of "geometry" to get flux into box
  EkIn=EkIn.*r*geometry(sect);

  % add to EkmanIn
  for stn=1:length(lon2);
    dg=(gamma_ek(stn)-glevels).^2;
    dg=find(dg==min(dg));
    EkmanInM(dg,sect)=EkmanInM(dg,sect)+EkIn(stn);
    EkmanInH(dg,sect)=EkmanInH(dg,sect)+EkIn(stn)*temp_ek(stn);
    EkmanInS(dg,sect)=EkmanInS(dg,sect)+EkIn(stn)*salam_ek(stn);
  end;

end;	% loop over sections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%