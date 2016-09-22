function prepctd(filedir,filename,ctdsal,ctdptemp,ctdpres,bin_press,lat,lon,dens_choice,file_type,meth_bot_tri) 
%==========================================================================
% PREPCTD     Prepares the CTD data to be incorporate into DIABOX with 
%             appropriate variable names and dimensions.
%
% The initial data processing is user-specific given how the initial ctd
% data format will vary according to the choice of processing done during
% the cruise. To simplify, here are the requirements to set up the initial
% variables that will be used by the inverse model.
% DIABOX requires temperature, salinity and pressure on a regular grid
% (ndepths,nstations):
%   In situ temperature, salinity and pressure matrices [ndepths,nstations]
%   Latitude and longitude of the stations [1,nstations]
%
% DESCRIPTION
% This function performs a number of tasks:
%     1) Check consistency of the input data
%     2) Calculate the area of the bottom triangles
%     3) Calculate the pressure of the isopycnal surfaces chosen. So far,
%     only two choices: Neutral density and potential density. This part of
%     the code needs to be updated according to the user's choice of
%     isopycnal surfaces.
%     4) Calculates geostrophic velocity referenced to the surface between  
%        stations to deepest common depth.
%     5) Adds the geostrophic velocity to area beneath deepest 
%        common depth. 
%     6) Extrapolate temperature and salinity down into the bottom triangle
%     7) Calculate salinity anomalies
%
% INPUTS:
%    - filedir: Path to the CTD data
%    - filename: CTD section name.
%    - ctdsal:  Name of salinity from CTD cast    [n*m]
%    - ctdptemp: Name of pot. temperature from CTD cast [n*m]
%    - ctdpres: Name of pressure from CTD cast    [n*m]
%    - bin_press: Pressure grid                   [n*1]
%    - lat: Name of latitude                      [1*m]
%    - lon: Name of longitude                     [1*m]
%    - dens_choice: Density to be used (neutral 
%                   density or pot. den           [str]
%    - meth_bot_tri: Method choice to extrapolate the geostrophic velocity 
%    in the bottom triangles. 2 choices: 
%           1) 'pad': Latest good value is padded down to the bottom. 
%           2) 'phi': extrapolation of the geopotential anomaly shear 
%                    following Thompson and Heywood, 2008. 
%
% OUTPUTS: 3 Matlab files are stored in "sect_dir":
%    - filename_raw: bin_ptemp bin_sal bin_press surf_press binPair_vel Ln 
%                    lat lon bottom frac
%    - filename_ns: glevels gtype
%    - filename_tria: area_tri
% 
%------------
% prepctd.m - last modified Apr. 2016
%
%==========================================================================
%% Load variables
if strcmp(file_type,'.nc') == 1
    disp('Opening netcdf files')
    eval(['ncid = netcdf.open(''', filedir , filename , file_type ,''');'])
    eval(['varID = netcdf.inqVarID(ncid,''' ctdsal ''');']);
    bin_sal = netcdf.getVar(ncid,varID,'double');
    eval(['varID = netcdf.inqVarID(ncid,''' ctdptemp ''');']);
    bin_ptemp = netcdf.getVar(ncid,varID,'double');
    eval(['varID = netcdf.inqVarID(ncid,''' ctdpres ''');']);
    pressmat = netcdf.getVar(ncid,varID,'double');
    eval(['varID = netcdf.inqVarID(ncid,''' bin_press ''');']);
    bin_press = netcdf.getVar(ncid,varID,'double');
    eval(['varID = netcdf.inqVarID(ncid,''' lat ''');']);
    lat = netcdf.getVar(ncid,varID,'double');
    eval(['varID = netcdf.inqVarID(ncid,''' lon ''');']);
    lon = netcdf.getVar(ncid,varID,'double');
    netcdf.close(ncid)
elseif strcmp(file_type,'.mat') == 1
    disp('Opening matlab files')
    eval(['load ',  filedir , filename,file_type,' ',ctdsal,' ',ctdtemp,' ',ctdpres, ' ',lat,' ',lon,'']);
else
    error('Please choose either a .nc or .mat')
end
%% Initialise variables and check data are consistent
disp('processing data')
% Check dimensions.
if size(lat,2) == 1; lat = lat'; end;
if size(lon,2) ==1; lon = lon'; end;

np = size(ctdpres,1); % Size of the ctd pressure grid

nstations = length(lat); % Number of stations
nstatPair = nstations - 1; % Number of station pairs.

binPair_vel = ones(np,nstatPair).*NaN;

%...Check that salinity, temperature and pressure
%...values are consistent.

fs=sum(~isnan(bin_sal));
ft=sum(~isnan(bin_ptemp));
ms=fs-ft;

ind = find(ms~=0);
if isempty(ind)~=1
 disp('Still some inconsistencies in the data')
 disp('This should not be the case since date have been interpolated');
 disp('Check the data again');
 return
end
%... Check the reference frame of longitude. 
%... For SW_prop and dobox longitude -180...+180 but to 
%... calculate neutral surface information lon =0...360.
%... Latitude and longitude as row vector.
[mlat,~]= size(lat);
 if mlat == 1
     lat=lat;
 else 
     lat=lat';
     lon=lon';
 end
xc=find(lon > 180);
lon(xc)=lon(xc)-360;

%% Calculate area of bottom triangle 
[~,bot_bin] = lastgood(bin_sal);
for istat = 1:nstations
    bottom(istat) = pressmat(bot_bin(istat),istat);
end
max_dpth = sw_dpth(bottom, lat);

%% Calculate pressure of density surfaces
if strcmp(dens_choice,'gamma') % Neutral density
    glev;
    sns=NaN*ones(length(glevels),length(lon));
    tns=sns;pns=sns;dsns=sns;dtns=sns;dpns=sns;
    disp('Calculating pressure of neutral surfaces.');
    for i=1:length(lon);
        fprintf('  Station %d of %d.\n',i,length(lon));
        [snst,tnst,pnst,dsnst,dtnst,dpnst] = ...
            neutral_surfaces(bin_sal(:,i),bin_temp(:,i), ...
        pressmat(:,i),gamma(:,i),glevels);
        sns(:,i)=snst; tns(:,i)=tnst; pns(:,i)=pnst;
        dsns(:,i)=dsnst; dtns(:,i)=dtnst; dpns(:,i)=dpnst;
    end;
    clear snst tnst pnst dsnst dtnst dpnst;
    sns=change(sns,'==',-99.0,NaN);
    tns=change(tns,'==',-99.0,NaN);
    pns=change(pns,'==',-99.0,NaN);
    dsns=change(dsns,'==',-99.0,NaN);
    dtns=change(dtns,'==',-99.0,NaN);
    dpns=change(dpns,'==',-99.0,NaN);

    surf_press=pns;
    
    gtype = []; % Type of density used. if empty = Neutral density.
    save glevels.mat glevels gtype;
    
elseif strcmp(dens_choice,'pden') % Potential density
    
    % This part of the code needs to be updated according to which
    % isopycnals will be used for the inversion.
    % calculate density referenced to different depths
    bin_temp = sw_temp(bin_sal,bin_ptemp,bin_press,0);
    
    sig0    = sw_dens0(bin_sal,bin_temp)-1000.;
    sig500  = sw_pden(bin_sal,bin_temp,bin_press,500)-1000.;
    sig1000 = sw_pden(bin_sal,bin_temp,bin_press,1000)-1000.;
    sig1500 = sw_pden(bin_sal,bin_temp,bin_press,1500)-1000.;
    sig2000 = sw_pden(bin_sal,bin_temp,bin_press,2000)-1000.;
    sig2500 = sw_pden(bin_sal,bin_temp,bin_press,2500)-1000.;

    disp('Calculating pressure of potential density surfaces');
    
    layer_s0  = [ 28.9 29.06 29.11 ];
    layer_s500 = [  ];
    layer_s1000 = [];
    layer_s1500 = [  ];
    layer_s2000 = [ ];
    layer_s2500 = [ ];
    
    glevels = [layer_s0'];
    gtype   = [0;0;0]; % Type of potential density used. 0 = referenced 
                          % to the surface. 1000 = referenced to 1000m
                        
    save glevels.mat glevels gtype;
    
    %layer_s25 = [ 39.744 39.746 ];
    if isempty(layer_s0) ~=1
       [surf_prs0,~,~,~]  = calwmb(layer_s0,sig0,pressmat);
    else
        surf_prs0 = [];
    end
    
    if isempty(layer_s500) ~=1
       [surf_prs500,~,~,~] = calwmb(layer_s500,sig500,pressmat);
    else
       surf_prs500 = [];
    end
    
    if isempty(layer_s1000) ~=1
       [surf_prs1000,~,~,~] = calwmb(layer_s1000,sig1000,pressmat);
    else
       surf_prs1000 = 0;
    end
            
    if isempty(layer_s1500) ~=1
       [surf_prs1500,~,~,~] = calwmb(layer_s1500,sig1500,pressmat);
    else
       surf_prs1500 = [];
    end
    
    if isempty(layer_s2000) ~=1
       [surf_prs2000,~,~,~] = calwmb(layer_s2000,sig2000,pressmat);
    else
       surf_prs2000 = [];
    end
    
    if isempty(layer_s2500) ~=1
       [surf_prs2500,~,~,~] = calwmb(layer_s2500,sig2500,pressmat);
    else
       surf_prs2500 = [];
    end

    % summarize the layer_press.
    surf_press = [surf_prs0];

else
    error('Please choose the correct density choice (''gamma'' or ''pden'')');
end
% % added 28 May 2002: for each cast, assign
% % deepest existing level+1 to bottom depth and
% % shallowest existing level-1 to surface.
% for i=1:length(lon);
%   q=find(isfinite(surf_press(:,i)));
%   if (min(q)>1)
%     surf_press(min(q)-1,i)=0;
%   end;
%   if (max(q)<length(glevels));
%     surf_press(max(q)+1,i)=bin_press(index(i));
%   end;
%   if ~length(q);	% only one layer present
% 			% (and thus no levels)
%     q2=min(find(glevels>=nanmean(gamma(:,i))));
%     surf_press(q2,i)=bottom(i);
%     if (q2>1)
%       surf_press(q2-1,i)=0;
%     end;
%   end;
% end;

%% Geostrophic velocities
%...DESCRIPTION
%...This file uses Phil Morgan's Seawater routines (sw_gvel, sw_gpan,
%...sw_dpth, and sw_dist). It calculates the geostrophic velocity
%...relative to the sea surface (sw_gvel). The velocity is then calculated
%...relative to the deepest common depth between each cast.
bin_temp = sw_temp(bin_sal,bin_ptemp,pressmat,0);
gpan     = sw_gpan(bin_sal,bin_temp,pressmat);
gvel     = sw_gvel(gpan,lat,lon);

%% Calculate the bottom triangles. 
% Two options:
% 'pad' - The simplest: Just extrapolate the velocity at the DCL down to
% the bottom.
% 'phi' - Calculate the geostrophic velocities in the bottom triangle by
% extrapolating the geopotential anomaly as in Thompson and Heywood, 2008.

if strcmp(meth_bot_tri,'pad');
   % just pad gvel to deepest level
   gvel2 = gvel;
   [gvelp,frac,area_tri] = gvelpad(bottom,gvel,gvel2,bin_press,lat,lon);
elseif strcmp(meth_bot_tri,'phi');
   gvel2 = sw_bottri_gpan(gpan,lon,lat);
   [gvelp,frac,area_tri] = gvelpad(bottom,gvel,gvel2,bin_press,lat,lon);
else
  disp('invalid entry for extrapolation method - gvel not calculated')
end

%% pad temp and salinity to corrpress
[salp]   = padit(bin_sal,pressmat,bottom); 
[ptempp] = padit(bin_ptemp,pressmat,bottom); 

%% Create variables
binPair_vel = gvelp;
bin_sal     = salp;
bin_potmp   = ptempp;

% Salinity anomaly  (!!!This can be changed but it also has to be changed
% in prepctd.m!!!)
bin_salam = (bin_sal - 35)*1e-3;

% calculate 2D field of layer numbers
Ln=laynumfield(surf_press,ctdpres,bottom);

%% Save appropriate variable to .mat files
load dir_loc.mat sect_dir
savearray1= ' bin_potmp bin_sal bin_salam bin_press surf_press binPair_vel Ln lat lon bottom frac';
savearray3= ' glevels gtype';
savearray4= ' area_tri';

eval(['save ', sect_dir filename '_raw', savearray1,';']);
eval(['save ', sect_dir filename '_ns', savearray3,';']);
eval(['save ', sect_dir filename '_tria', savearray4,';']);

disp('Files saved - prepctd finished')

return
%% Subroutines
function [surf_prs,avewmb,stdwmb,nstnwmb] = calwmb(defsig,sig,pressmat)

% CALWMB  Calculate water mass boundary
%=========================================================================
%
% USAGE:  surf_prs = calwmb(defsig,sig,bin_press)
%
% DESCRIPTION:
%    Calculate pressure which corresponds to defined sigma surface
%
% INPUT:  (all must have same dimensions)
%   defsig = sigma surface setted to boundary [kg/m^3 ex. sig0=25.4]
%
%   sig    = sigma field; [ndep,nstn]=size(sig) [kg/m^3]
%   bin_press = pressure; [1,ndep]=size(bin_press)
%
% OUTPUT:
%   surf_prs = pressure of water mass boundary  [pressure]
%
%=========================================================================

% ------ example of input --------
% load fram2005_prep.mat
% ptmp=sw_ptmp(bin_sal,bin_temp,bin_press,0.);
% % for sig0
% defsig = [ 27.7 27.97 ];
% sig=sw_dens(bin_sal,ptmp,0.)-1000.;

% for sig15
%defsig = [ 35.126 35.142 ];
%sig=sw_dens(bin_sal,ptmp,1500)-1000.;
%-----------------------------------


ndefsig  = length(defsig); % calc number of boundary
[ndep,nstn]=size(sig);

% prepare the matrix
surf_prs=ones(ndefsig,nstn).*0.*inf;

% Calculate pressure which corresponds to defined sigma surface

for Jwmb=1:ndefsig; % loop for water mass boundary
    threpr = defsig(Jwmb);
    for Istn=1:nstn; % loop for each station
        if( any (sig(:,Istn) >= threpr) == 1 );
            [jj,ii] = find ( sig(:,Istn) >= threpr );
            surf_prs(Jwmb,Istn)=pressmat(jj(1),Istn);
        end;
    end;
end

% Calculate average density layer

for JJ=1:ndefsig; % loop for each layer
    A = surf_prs(JJ,:);
    A = change(A,'==',0,NaN  ); % I am not sure this is neccesary or not for general purpose.
    avewmb(JJ)=mean(A(~isnan(A)));
    stdwmb(JJ)=std(A(~isnan(A)));
    nstnwmb(JJ)=length(A(~isnan(A)));
end


 %%
% figure;
% hold on
% 
% % temp sec
% pcolor(lon,bin_press,bin_temp); axis ij; shading flat;
% colorbar; caxis([-4 6]);
% 
% % Defined WM boundary
% plot(lon,surf_prs','k','LineWidth',3);
% 
% hold off
% axis ij;
% title('temp sect with surfpress(dot)');
% xlabel('longitude');ylabel('pressure')

return

