%--------------------------------------------------------------------------
% Read the Dai and Trenberth Global River Flow and Continental Discharge 
% Dataset
%
% Refs.: Dai,  A., and K. E. Trenberth,  2002: Estimates of freshwater 
% discharge from continents: Latitudinal and seasonal variations. 
% J. Hydrometeorol., 3, 660-687
%--------------------------------------------------------------------------
clear all

lat_river  = nc_varget('coastal-stns-Vol-monthly.updated-oct2007.nc','lat_mou');
lon_river  = nc_varget('coastal-stns-Vol-monthly.updated-oct2007.nc','lon_mou');
river      = nc_varget('coastal-stns-Vol-monthly.updated-oct2007.nc','riv_name');
river_flow = nc_varget('coastal-stns-Vol-monthly.updated-oct2007.nc','FLOW');


for iriver = 1:length(lat_river)
    R(iriver) = nanmean(river_flow(:,iriver)); 
end

q = find(isnan(R));

R(q) = [];
lat_river(q) = [];
lon_river(q) = [];
river(q,:)  = [];

% add silica flux (kmol/s)
si_conc=9.6*ones(size(R));      % SiO_2 (mg/L)
for i=1:length(R);
 if (strcmp(char(river(i)),'Amazon'))
   si_conc(i)=7.2;
 elseif (strcmp(char(river(i)),'Orinoco'))
   si_conc(i)=11.5;
 elseif (strcmp(char(river(i)),'Zaire'))
   si_conc(i)=9.8;
 elseif (strcmp(char(river(i)),'Mississippi'))
   si_conc(i)=7.6;
 end;
end;

si_flux=si_conc.*(R*100/1e6)/6.00843;

save river_runoff.mat lat_river lon_river R si_flux river