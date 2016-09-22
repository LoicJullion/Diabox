function [lon_era,lat_era,lonvec,latvec,sss,index,sst,netheat,latent,sensible,longwave,shortwave,eminusp,evap,precip,taux,tauy]=loaderayear(year)
%==========================================================================
% function [lon_era,lat_era,lonvec,latvec,sss,index, ...
%        sst,netheat,eminusp,taux,tauy]=loaderayear(year,clim_asdir,climdir);
%
% load a year (=1979 through 2012) of monthly observations from the 
% ERA-Interim reanalysis fields.
%
% INPUT:
%   - year: The year(s) you want to use. -99 = climatology
%   - clim_asdir: Path to the ERA-Interim reanalysis
%   - climdir: Path to the oceanic climatology. We use MIMOC 
%
% AUTHORS: Loic Jullion based on Rick Lumpkin's code.
% 
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
load dir_loc.mat climdir clim_asdir 
    %% Load the data
    
    fileext='.cdf';
    
    eval(['ncid = netcdf.open(''' clim_asdir 'sst_erai.cdf'');']);
    varID   = netcdf.inqVarID(ncid,'longitude');
    lon_era = netcdf.getVar(ncid,varID);
    
    % break data at dateline instead of prime meridian
    q=find(lon_era>=180);
    lon_era(q)=lon_era(q)-360;
    [lon_era,ind_lon]=sort(lon_era);
    
    varID   = netcdf.inqVarID(ncid,'latitude');
    lat_era = netcdf.getVar(ncid,varID);
    netcdf.close(ncid);
    
    if (year<0)
        disp('Loading ERA-Interim long-term monthly means');
        eval([' load ' clim_asdir 'era_interim_monthly_clim.mat; ']);
        netheat   = netheat_mn_clim;
        precip    = precip_mn_clim;
        evap      = evap_mn_clim;
        sst       = sst_mn_clim;  
        taux      = taux_mn_clim;
        tauy      = tauy_mn_clim;
        longwave  = longwave_mn_clim;
        shortwave = shortwave_mn_clim;
        latent    = latent_mn_clim;
        sensible  = sensible_mn_clim;

    else
      
      fprintf('Loading ERA-Interim fields, year %d.\n',year);
      eval(['ncid = netcdf.open(''' clim_asdir 'sst_erai.cdf'');']);
      varID = netcdf.inqVarID(ncid,'time');
      t     = netcdf.getVar(ncid,varID,'double');
      netcdf.close(ncid);
      
      % time: hours since 1900/1/1 00:00:00
      % For the monthly means of daily means, the time is: [Year,Month,Day,0,0,0]
      t = t.*(60*60);
      q = find(t>=date2unixsecs(year,1,1,0,0,0,1900) & ...
                t<=date2unixsecs(year,12,31,0,0,0,1900));
      corner = [0 0 q(1)-1];
      endpt  = [length(lon_era) length(lat_era) length(q)];

      % For synoptic monthly means, the time is different: The first synoptic
      % mean is at time [Year,Month,Day,12,0,0] and the second synoptic mean is
      % at time [Year,Month,Day+1,0,0,0]
            
      eval(['ncid = netcdf.open(''' clim_asdir 'latent_hflx_eraI' fileext ''');']);
      varID = netcdf.inqVarID(ncid,'time'); 
      time2 = netcdf.getVar(ncid,varID,'double');
      netcdf.close(ncid);
  
      time2=time2*(60*60);
      q2=find(time2>=date2unixsecs(year,1,1,0,0,0,1900) & ...
                time2<=date2unixsecs(year,12,31,0,0,0,1900));

      corner_fcst = [0 0 q2(1)-1];
      endpt_fcst  = [length(lon_era) length(lat_era) length(q2)];
   
      % load, changing all heat flux signs to retain
      % convention that POSITIVE = heat INTO ocean
      % Montlhy synoptic means of accumulated forecast        
      eval(['ncid = netcdf.open(''' clim_asdir 'latent_hflx_eraI' fileext ''');']);
      varID   = netcdf.inqVarID(ncid,'slhf'); 
      latent2 = netcdf.getVar(ncid,varID,corner_fcst,endpt_fcst,'double');
      netcdf.close(ncid);
      
      eval(['ncid = netcdf.open(''' clim_asdir 'sensible_hflx_erai' fileext ''');']);
      varID   = netcdf.inqVarID(ncid,'sshf'); 
      sensible2 = netcdf.getVar(ncid,varID,corner_fcst,endpt_fcst,'double');
      netcdf.close(ncid);
  
      eval(['ncid = netcdf.open(''' clim_asdir 'longwave_hflx_erai' fileext ''');']);
      varID   = netcdf.inqVarID(ncid,'str'); 
      longwave2 = netcdf.getVar(ncid,varID,corner_fcst,endpt_fcst,'double');
      netcdf.close(ncid);

      eval(['ncid = netcdf.open(''' clim_asdir 'shortwave_hflx_erai' fileext ''');']);
      varID   = netcdf.inqVarID(ncid,'ssr'); 
      shortwave2 = netcdf.getVar(ncid,varID,corner_fcst,endpt_fcst,'double');
      netcdf.close(ncid);

      eval(['ncid = netcdf.open(''' clim_asdir 'tot_precipitation_erai' fileext ''');']);
      varID   = netcdf.inqVarID(ncid,'tp'); 
      precip = netcdf.getVar(ncid,varID,corner_fcst,endpt_fcst,'double');
      netcdf.close(ncid);

      eval(['ncid = netcdf.open(''' clim_asdir 'evap_erai' fileext ''');']);
      varID = netcdf.inqVarID(ncid,'e'); 
      evap  = netcdf.getVar(ncid,varID,corner_fcst,endpt_fcst,'double');
      netcdf.close(ncid);


      % Montlhy means of daily means    
      eval(['ncid = netcdf.open(''' clim_asdir 'sst_erai' fileext ''');']);
      varID = netcdf.inqVarID(ncid,'sst'); 
      sst   = netcdf.getVar(ncid,varID,corner,endpt,'double');
      netcdf.close(ncid);
      % convert from K to degC
      sst=sst-273;

      
      eval(['ncid = netcdf.open(''' clim_asdir 'u_stress_erai' fileext ''');']);
      varID = netcdf.inqVarID(ncid,'iews'); 
      taux   = netcdf.getVar(ncid,varID,corner,endpt,'double');
      netcdf.close(ncid);

      eval(['ncid = netcdf.open(''' clim_asdir 'v_stress_erai' fileext ''');']);
      varID  = netcdf.inqVarID(ncid,'inss'); 
      tauy   = netcdf.getVar(ncid,varID,corner,endpt,'double');
      netcdf.close(ncid);

      % define netheat
      netheat=latent2+sensible2+longwave2+shortwave2;

      clear latent sensible longwave shortwave;
  
      latent    = latent2;
      longwave  = longwave2;
      shortwave = shortwave2;
      sensible  = sensible2;
    
      % re-order variables so longitude is -180:180
      netheat   = netheat(ind_lon,:,:);
      sst       = sst(ind_lon,:,:);
      taux      = taux(ind_lon,:,:);
      tauy      = tauy(ind_lon,:,:);
      latent    = latent(ind_lon,:,:);
      longwave  = longwave(ind_lon,:,:);
      shortwave = shortwave(ind_lon,:,:);
      sensible  = sensible(ind_lon,:,:);
      evap      = evap(ind_lon,:,:);
      precip    = precip(ind_lon,:,:);
    
    end
    %% calculate eminusp:
    % precip is always positive
    % evap is mostly negative
    % eminusp > 0 => evaporation dominates, 
    % eminusp < 0 => precip dominates.
    eminusp=-evap-precip;

    %%
    % reshape data to have the same format as
    % the COADS fields (time,index)
    
    eval(['ncid = netcdf.open(''' clim_asdir 'Mask' fileext ''');']);
    varID  = netcdf.inqVarID(ncid,'lsm'); 
    landmask   = netcdf.getVar(ncid,varID,'double');
    netcdf.close(ncid);
    landmask = landmask(ind_lon,:);
    clear ind_lon;

    % redefine so landmask=1 over ocean, NaN over land
    landmask = landmask+1;
    % Replace grid points on land by NaNs,
    landmask(find(landmask==2))=NaN;
    % Get indices of all the grid points
    index=1:size(landmask,1)*size(landmask,2);
    % Transform the landmask in a vector
    landmask=reshape(landmask, ...
        [1 size(landmask,1)*size(landmask,2)]);

    % Create matrices containing lon and lat
    lonvec=lon_era*ones(1,length(lat_era)); % Lon matrix
    latvec=ones(length(lon_era),1)*lat_era';  % Lat matrix
    
    mat   = NaN.*landmask;
    
    sst2       = repmat(mat,12,1); netheat2  = repmat(mat,12,1);
    taux2      = repmat(mat,12,1); tauy2     = repmat(mat,12,1);
    latent2    = repmat(mat,12,1); sensible2 = repmat(mat,12,1);
    shortwave2 = repmat(mat,12,1); longwave2 = repmat(mat,12,1);
    eminusp2   = repmat(mat,12,1); precip2   = repmat(mat,12,1);
    evap2      = repmat(mat,12,1);
    
    % Vectorize
    lonvec=reshape(lonvec,[size(mat,1) size(mat,2)]);
    latvec=reshape(latvec,[size(mat,1) size(mat,2)]);

    for imonth =1:12;    
    netheat2(imonth,:) = reshape(squeeze(netheat(:,:,imonth)),[size(mat,1) ...
                                 size(mat,2)]);

    longwave2(imonth,:) = reshape(squeeze(longwave(:,:,imonth)),[size(mat,1) ...
                                 size(mat,2)]);
                             
    shortwave2(imonth,:) = reshape(squeeze(shortwave(:,:,imonth)),[size(mat,1) ...
                                 size(mat,2)]);
                             
    latent2(imonth,:) = reshape(squeeze(latent(:,:,imonth)),[size(mat,1) ...
                                 size(mat,2)]);
                             
    sensible2(imonth,:) = reshape(squeeze(sensible(:,:,imonth)),[size(mat,1) ...
                                 size(mat,2)]);
                            
    sst2(imonth,:)     = reshape(squeeze(sst(:,:,imonth)),[size(mat,1) ...
                                 size(mat,2)]);

    taux2(imonth,:)    = reshape(squeeze(taux(:,:,imonth)),[size(mat,1) ...
                                 size(mat,2)]);

    tauy2(imonth,:)    = reshape(squeeze(tauy(:,:,imonth)),[size(mat,1) ...
                                 size(mat,2)]);

    eminusp2(imonth,:) = reshape(squeeze(eminusp(:,:,imonth)),[size(mat,1) ...
                                 size(mat,2)]);
    precip2(imonth,:) = reshape(squeeze(precip(:,:,imonth)),[size(mat,1) ...
                                 size(mat,2)]);
    evap2(imonth,:) = reshape(squeeze(evap(:,:,imonth)),[size(mat,1) ...
                                 size(mat,2)]);

    end
    
    clear sst eminusp netheat taux tauy longwave shortwave latent sensible evap precip
    
    sst       = sst2; 
    eminusp   = eminusp2;
    taux      = taux2;
    tauy      = tauy2;
    netheat   = netheat2;
    latent    = latent2;
    longwave  = longwave2;
    shortwave = shortwave2;
    sensible  = sensible2;
    precip    = precip2;
    evap      = evap2;

    % retain values over the ocean
    q       = find(isfinite(landmask));
    index   = index(q);
    lonvec  = lonvec(q);
    latvec  = latvec(q);
    netheat = netheat(:,q);
    latent = latent(:,q);
    sensible = sensible(:,q);
    longwave = longwave(:,q);
    shortwave = shortwave(:,q);
    sst     = sst(:,q);
    taux    = taux(:,q);
    tauy    = tauy(:,q);
    eminusp = eminusp(:,q);
    evap    = evap(:,q);
    precip  = precip(:,q);

    clear landmask

% load MIMOC Sea surface salinity data.
eval(['load ',climdir,'sss.mat sss; ']);

return

