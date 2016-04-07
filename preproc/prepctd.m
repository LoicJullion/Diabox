function prepctd(filename,ctdsal,ctdtemp,ctdpres,gamma,lat,lon) 
%==========================================================================
% PREPCTD     Prepares the CTD data to be incorporate into DIABOX with 
%             appropriate variable names and dimensions.
%
% Ideally, DIABOX likes data interpolated on a regular grid from 0 to the
% bottom with no duplicate depths and pressure inversions.
% The best option is probably to check the data before running them
% through this function to make sure the data are fine.
%
% DESCRIPTION
% This file does a number of tasks 
%     1) Test to see that a temperature and salinity values exist
%        at every 2db data point within a cast. If not 
%        it interpolates to find these values.
%     2) Calculates geostrophic velocity between stations to 
%        deepest common depth.
%     3) Adds the geostrophic velocity to area beneath deepest 
%        common depth using (v=const. ie. dvdz=0). Calculates 
%        the area of the bottom triangle 
%     4) Find temperature and salinity values at standard
%        depths. In the upper 500 db this is the actual value
%        property value at the standard depth. Below this it is 
%        a 20 db mean centred around standard depth.
%     5) Calculates the gamma_n values for the section, using the 
%        2db data. Then finds the temperature, salinity and 
%        pressure value for chosen neutral surfaces. The mean 
%        pressure between these stations of neutral surfaces 
%        is then calculated.
%
% INPUTS:
%    - filename: section name
%    - ctdsal:  Name of salinity from CTD cast    [n*m]
%    - ctdtemp: Name of temperature from CTD cast [n*m]
%    - ctdpres: Name of pressure from CTD cast    [n*m]
%    - gamma: Name of density from CTD cast       [n*m]
%    - lat: Name of latitude                      [1*m]
%    - lon: Name of longitude                     [1*m]
%
%
% OUTPUTS:
% 
%------------
% prepctd.m - last modified Mar. 2014
%
%==========================================================================
disp('processing data')
% Check dimensions.
if size(lat,2) == 1; lat = lat'; end;
if size(lon,2) ==1; lon = lon'; end;
if size(ctdpres,1) == 1; ctdpres = ctdpres'; end
%... check that ctdpres is monotonic

% Decide on the size (in the vertical) you want all you sections to be
% equal to
vert_size = 7000; % in db

% Standard pressure grid
dp = 2; % in db
% Pressure grid from 0 to vert_size every dp
bin_press =[0:dp:vert_size]';

d = size(ctdpres,1); % Size of the ctd pressure grid
nstations = length(lat); % Number of stations
nstatPair = nstations - 1; % Number of station pairs.
% Declare variables
bin_temp = ones(length(bin_press),nstations).*NaN;
bin_sal  = bin_temp;

binPair_vel = ones(length(bin_press),nstatPair).*NaN;

% For each profile, loop through casts.
% Check consistency between pressure, temp and sal.
% Make sure there are no pressure inversion or duplicate pressure.
% Interpolate the ctd data on the bin_press grid to create bin_sal and
% bin_temp
for i=1:nstations;
    pp = ctdpres(:,i);
    tt = ctdtemp(:,i);
    ss = ctdsal(:,i);
    
    gd_pp = find(isnan(pp)~=1);
    gd_tt = find(isnan(tt)~=1);
    gd_ss = find(isnan(ss)~=1);
    
    % make sure all the data have data at the same place.
    if (length(gd_pp) ~= length(gd_ss))
        disp('Not the same number of NaNs in the profiles of press and sal');
    elseif length(gd_pp) ~= length(gd_tt)
        disp('Not the same number of NaNs in the profiles of press and temp');
    elseif length(gd_ss) ~= length(gd_tt)
        disp('Not the same number of NaNs in the profiles of temp and sal');
    end
    % pressure increment between each bin
    res_p = diff(pp);
    
    ind1 = find(res_p>0); % consecutive bins with increasing pressure    
    ind2 = find(res_p==0); % consecutive bins at the same pressure
    ind3 = find(res_p<0); % consecutive bins with reversed pressure
    
    if isempty(ind3) ~= 1
       fprintf(1,'Station %i: \n',i); 
       disp('Pressure is not monotonically increasing. Sorting the profile in descending order');
       [pp,pos] = sort(pp);
       tt       = tt(pos);
       ss       = ss(pos);
       clear pos;
    end    
    
    if isempty(ind2) ~=1
       fprintf(1,'Station %i: \n',i); 
       fprintf(1,'Bins with same pressure, get rid of the extra values\n');
       [~,Ia,~] = unique(pp(gd_pp));
       pp = pp(Ia);
       tt = tt(Ia);
       ss = ss(Ia);
    end

    % Refind the indices with no NaNs:
    gd_pp = find(isnan(pp)~=1);
    gd_tt = find(isnan(tt)~=1);
    gd_ss = find(isnan(ss)~=1);
    
    % Now that we are sure all the data are sorted in ascending order and
    % that all bins with the same pressure have been removed, interpolate
    % on the bin_press grid.
    bin_temp(:,i) = interp1(pp(gd_pp==gd_tt),tt(gd_pp==gd_tt),bin_press,'nearest');
    bin_sal(:,i)  = interp1(pp(gd_pp==gd_ss),ss(gd_pp==gd_ss),bin_press,'nearest');
    % Find the position of the first and last good values.
    [~,top_ind] = firstgood(bin_temp(:,i));
    [~,bot_ind] = lastgood(bin_temp(:,i));
    clear ind;
    ind = find(isnan(bin_sal(top_ind:bot_ind,i))==1);
    % In theory, since the data have been interpolated, there should be no
    % NaNs.
    if isempty(ind)~=1
        fprintf(1,'NaNs are present in the profile. Check the data\n');
    end
    clear ind;
    if top_ind~=1
        disp('Missing data near the surface');
        disp('Average the first 4 data points');
        bin_temp(1:top_ind-1,i) = mean(bin_temp(top_ind:top_ind+3,i));
        bin_sal(1:top_ind-1,i) = mean(bin_sal(top_ind:top_ind+3,i));
    end
end
% Now binned properties have been calculated

%...Check that salinity, temperature and pressure
%...values are consistent.

fs=sum(~isnan(bin_sal));
ft=sum(~isnan(bin_temp));
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


%------------------------
% June 2015: Not sure the following lines are still needed.
% I'll leave them just in case but they can probably be removed
% fix gaps in bin_temp which will produce
% "binPair_vel not consistent with bin_salam"
% in DIABOX
q1=find(~isnan(bin_temp));
q2=find(~isnan(bin_sal));
if (length(q1)<length(q2))
  q=find(~isnan(bin_sal) & isnan(bin_temp));
  for i=1:length(q);
    a=bin_temp(q(i)-1)-bin_temp(q(i)-2);
    b=bin_temp(q(i)-1)-2*a;
    bin_temp(q(i))=3*a+b;
  end
end
%------------------------ 
%...Calculate area of bottom triangle with respect to standard depths.
[~,bot_bin] = lastgood(bin_sal);
bottom = bin_press(bot_bin);
max_dpth = sw_dpth(bottom', lat);

dy=0.5 *abs(diff(max_dpth));
dx=sw_dist(lat,lon,'km')*1000;
area_tri=(dx.*dy)*0.5;

pres = bin_press*ones(1,nstations);

glev;

sns=NaN*ones(length(glevels),length(lon));
tns=sns;pns=sns;dsns=sns;dtns=sns;dpns=sns;
disp('Calculating pressure of neutral surfaces.');
for i=1:length(lon);
  fprintf('  Station %d of %d.\n',i,length(lon));
  [snst,tnst,pnst,dsnst,dtnst,dpnst] = ...
    	neutral_surfaces(bin_sal(:,i),bin_temp(:,i), ...
	pres(:,i),gamma(:,i),glevels);
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

% added 28 May 2002: for each cast, assign
% deepest existing level+1 to bottom depth and
% shallowest existing level-1 to surface.
for i=1:length(lon);
  q=find(isfinite(surf_press(:,i)));
  if (min(q)>1)
    surf_press(min(q)-1,i)=0;
  end;
  if (max(q)<length(glevels));
    surf_press(max(q)+1,i)=bin_press(index(i));
  end;
  if ~length(q);	% only one layer present
			% (and thus no levels)
    q2=min(find(glevels>=nanmean(gamma(:,i))));
    surf_press(q2,i)=bottom(i);
    if (q2>1)
      surf_press(q2-1,i)=0;
    end;
  end;
end;


%...DESCRIPTION
%...This file uses Phil Morgan's Seawater routines (sw_gvel, sw_gpan,
%...sw_dpth, and sw_dist). It calculates the geostrophic velocity
%...relative to the sea surface (sw_gvel). The velocity is then calculated
%...relative to the deepest common depth between each cast.
gpan = sw_gpan(bin_sal,bin_temp,bin_press);
gvel = sw_gvel(gpan,lat,lon);

% Calculate the bottom triangles. Two options:
% 'pad' - The simplest: Just extrapolate the velocity at the DCL down to
% the bottom.
% 'phi' - Calculate the geostrophic velocities in the bottom triangle by
% extrapolating the geopotential anomaly as in Thompson and Heywood, 2008.

meth_bot_tri = 'pad';

if strcmp(method1,'pad');
   % just pad gvel to deepest level
   gvel2 = gvel;
   [gvelp,frac,area_tri] = gvelpad(bottom,gvel,gvel2,bin_press,lat,lon);
elseif strcmp(method1,'phi');
   gvel2 = sw_bottri_gpan(gpan,lon,lat);
   [gvelp,frac,area_tri] = gvelpad(bottom,gvel,gvel2,bin_press,lat,lon);
else
  disp('invalid entry for extrapolation method - gvel not calculated')
end


%pad temp and salinity to corrpress
[salp]  = padit(bin_sal,bin_press,bottom); 
[tempp] = padit(bin_temp,bin_press,bottom); 

%calculate potemp
bin_ptmp = sw_ptmp(salp,tempp,bin_press,0);

binPair_vel = gvelp;
bin_sal     = salp;
bin_temp    = tempp;

bin_salam = (bin_sal - 35)*1e-3;
press=bin_press*ones(1,length(lon));
disp('Calculating potential temperature.');
bin_potmp = sw_ptmp(bin_sal,bin_temp,press,zeros(size(press)));
clear press;

% calculate sigma_theta
bin_sigthet=NaN*ones(size(bin_potmp));
bin_sigthet=sw_pden(bin_sal,bin_potmp, ...
	zeros(size(bin_sal)),zeros(size(bin_sal)))-1000;

% calculate 2D field of layer numbers
Ln=laynumfield(surf_press,ctdpres,bottom);

%... Prompts for file name and save appropriate variable to .mat files

savearray1= ' bin_temp bin_sal bin_press surf_press binPair_vel time';
savearray2= ' bin_salam bin_potmp bin_sigthet lat lon st_order bottom';
savearray22=' Ln directory';
savearray3= ' glevels gamma sns tns pns dsns dtns dpns';
savearray4= ' area_tri';

eval(['save ', filename '_raw', savearray1,savearray2,savearray22,';']);
eval(['save ', filename '_ns', savearray3,';']);
eval(['save ', filename '_tria', savearray4,';']);

disp('Files saved - prepctd finished')

