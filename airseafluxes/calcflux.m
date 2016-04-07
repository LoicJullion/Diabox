%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcflux.m
%
% calculate the density, salt and heat fluxes due to A/S forcing (from 
% fields of eminusp, netheat, sst, sss) and Ekman transport (from wind 
% stress fields taux, tauy).
%
% The air-sea fields are of dimension (month_num, index), with 
% corresponding latitude latvec and longitude lonvec. This routine
% calculates for each grid point of the vector the density, heat and salt
% fluxes due to A/S forcing.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% treat latitudes within +-lat_eq as a special
% region: f held constant from lat_eq to equator,
% and "Ekman" transport angle rotates from 90deg
% at lat_eq to 0deg (downwind) at equator
lat_eq=2;

%% Calculate gamma, alpha, beta and Cp 
% 1) Calculate gamma (the slow part!)
[nummonths,numdat]=size(sss);
disp('Calculating gamma.');
ssg=NaN*ones(nummonths,numdat);
latg=latvec;
q=find(latg>63.999);
if (length(q)) latg(q)=63.999; end;
q=find(latg<-79.999);
if (length(q)) latg(q)=-79.999; end;
long=lonvec;
q=find(long<0);
if (length(q)) long(q)=long(q)+360; end;
for month=1:nummonths;
  fprintf('  Month %d of %d.\n',month,nummonths);
  q=find(isfinite(sst(month,:)));
  if (length(q))
    ssg(month,:)=gamman(sss(month,:),sst(month,:),zeros(size(lonvec)),long,latg);
  end;
end;
q=find(ssg<-90);
ssg(q)=NaN;
clear long latg nummonths numdat

% 2) Calculate alpha (thermal expansion coeff., deg. C^-1), 
%    beta (saline contraction coeff., psu^-1) and Cp (specific heat, J/ kg C).
alpha = NaN*ones(size(sss));
beta  = alpha;
Cp    = alpha;
disp('Calculating alpha and beta.');
alpha = sw_alpha(sss,sst,zeros(size(sst)));
beta  = sw_beta(sss,sst,zeros(size(sst)));
disp('Calculating Cp.')
Cp    = sw_cp(sss,sst,zeros(size(sst)));

% calculate density ssd(0,S,0) (kg/m^3)
ssd   = NaN*ones(size(sss));
ssdS0 = ssd;
ssd   = sw_dens(sss,sst,zeros(size(sst)));
ssdS0 = sw_dens(zeros(size(sst)),sst,zeros(size(sst)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate density, heat and salt fluxes due to air-sea forcing.

% 1) Freshwater flux (m/s) (Positive = flux to ocean.
ffw          = NaN*ones(size(sss));
ffw          = -eminusp;
ffw_norivers = -orig_eminusp;

% error: error_multfact of E-P
ffwerr = error_multfact*abs(orig_eminusp);
clear orig_eminusp;

% 2) Heat flux (deg.C m/s) (Positive = flux to ocean)
fh = NaN*ones(size(sss));
fh = netheat./Cp./ssd;
% std. error of heat flux
%%%%%%%fherr=abs(netheaterr./Cp./ssd);
fherr = error_multfact*abs(fh);

% 3) Density flux (kg/m^2 s)
fm   = NaN*ones(size(sss));fmerr=fm;
fm   = ssd.*-alpha.*fh+ssdS0.*beta.*sss./(1-sss/1000).*eminusp;
fm_h = ssd.*-alpha.*fh;
% std. error of net density flux
fmerr = sqrt( (alpha.*ssd.*fherr).^2 + ...
              (ssdS0.*beta.*ffwerr.*sss./(1-sss/1000)).^2 );
clear alpha beta ssdS0 Cp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EKMAN COMPONENT:    
% calculate Aek (Sv), the net mass transport across density classes  
% caused by direct wind forcing		     
disp('Calculating Ekman flux');
f   = 2*7.2921e-05*sin(Latmat*2*pi/360);
Mek = NaN*ones(size(Latmat,1),size(Latmat,2), ...
	  size(sst,1));
Sek = Mek;Tek=Mek;
%%% TWEAK TO REMOVE SINGULARITY IN EKMAN TRANSPORT: f STAYS NONZERO AT EQUATOR
q = find(abs(Latmat)<lat_eq);
if (length(q))
    q2=find(abs(f)>0);
    f(q)=min(abs(f(q2)))*sign(Latmat(q)+1e-10);
end;

load glevels;
gamma2  = [glevels(1)-diff(glevels(1:2)); ...
           glevels;last(glevels)+last(diff(glevels))];
glayers = gamma2(1:length(gamma2)-1)+diff(gamma2)/2;

EkmanInM    = zeros(length(glayers),size(sectfiles,1));
EkmanInH    = EkmanInM; EkmanInS=EkmanInM;
DaysInMonth = [31 28.25 31 30 31 30 31 31 30 31 30 31];

% load Mixed layer climatological values of density, salinity,
% temperature and gamma from MIMOC, to calculate
% Ekman advection of properties
eval(['load ',climdir,'mld_mimoc_inversion.mat']);
q=find(lon_mimoc>=180);
lon_mimoc(q)=lon_mimoc(q)-360;
[lon_mimoc,ind_lon]=sort(lon_mimoc);

dens_mld = dens_mld(:,ind_lon,:);
temp_mld = temp_mld(:,ind_lon,:);
sal_mld = sal_mld(:,ind_lon,:);
gamm_mld = gamm_mld(:,ind_lon,:);
clear q ind_lon

for month=1:size(sst,1);
  fprintf('  month %d\n',month);
  Taux        = NaN*ones(size(Latmat));
  Tauy        = Taux;
  SST         = Taux;
  Rho         = Taux;
  Salam       = Taux;
  Temp        = Taux;
  Gamma       = Taux;
  Taux(index) = taux(month,:);
  Tauy(index) = tauy(month,:);
  SST(index)  = sst(month,:);

  i     = find(lat_mimoc>=min(lat_coads) & lat_mimoc<=max(lat_coads));
  j     = find(lon_mimoc>=min(lon_coads) & lon_mimoc<=max(lon_coads));
  Rho   = dens_mld(i,j,month-12*floor((month-1)/12));
  Salam = (sal_mld(i,j,month-12*floor((month-1)/12)) -35)/1e3;
  Temp  = temp_mld(i,j,month-12*floor((month-1)/12));
  Gamma = gamm_mld(i,j,month-12*floor((month-1)/12));
  
  % FLIP THE MATRIX. ISSUE WITH THE MIMOC DATA. THE VARIABLES NEED TO BE
  % SWITCHED UP AND DOWN, NOT SURE WHY. MIGHT BE DELETED IF THE ISSUE IS
  % RESOLVED OR ANOTHER DATA SET IS USED.
  Rho   = flipud(Rho);
  Salam = flipud(Salam);
  Temp  = flipud(Temp);
  Gamma = flipud(Gamma);
  % END FLIPPING
  
  % Ekman flux Ue and Ve (m^2/s)
  %%% fix to equatorial values where f -> 0:
  % Ekman flux is rotated an angle theta from
  % wind stress, theta rotates from 90 degrees
  % at 2N through 0 at equator to -90 at 2S.
  Ue       = -sqrt(-1)*(Taux+sqrt(-1)*Tauy)./Rho./f;
  magU     = abs(Ue);
  angU     = angle(Ue);
  angW     = angle(Taux+sqrt(-1)*Tauy);
  latmat   = lat_coads*ones(1,length(lon_coads));
  rot      = (lat_eq-abs(latmat)).*exp(sqrt(-1)*angW) + ...
              abs(latmat).*exp(sqrt(-1)*angU);
  ang      = angU;
  q        = find(abs(lat_coads)<=lat_eq);
  ang(q,:) = angle(rot(q,:));
  Ue       = magU.*cos(ang);
  Ve       = magU.*sin(ang);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % calculate Ekman fluxes into/out of box. %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Output has units of m^3/s times property units
  % (1 for volume, deg.C for temp, etc.)
  [EkM,EkH,EkS] = ekmanin(lon_coads,lat_coads, ...
        Ue,Ve,Gamma,Temp,Salam,glayers,file_geom,boxnumber);
  % add to previous months values, with a scaling
  % factor for the number of days in that month
  % (scalefact=1.0185 for 31d months, 0.9281 for February)
  scalefact = DaysInMonth(month-12*floor((month-1)/12))./mean(DaysInMonth);
  EkmanInM  = EkmanInM+EkM.*scalefact;
  EkmanInH  = EkmanInH+EkH.*scalefact;
  EkmanInS  = EkmanInS+EkS.*scalefact;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % calculate interior Ekman fluxes %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % retain values within box
  q       = find(isnan(SST));
  Rho(q)  = NaN; Salam(q) = NaN; 
  Temp(q) = NaN; Gamma(q) = NaN;
  Ue(q)   = NaN; Ve(q)    = NaN;

  % convert Ue, Ve from m^2/s to m^3/s
  [area,dx,dy] = sw_area(lon_coads,lat_coads);
  
  Ue = Ue.*dy;
  Ve = Ve.*dx;

  % extrapolate Ve, Ue to be one column/row
  % smaller than gridded field (so diff(Ue,1,2)
  % and diff(Ve,1,1) are on gridpoints).
  % make sure all edge values of extrapolation
  % are zero: cross-edge fluxes are in EkmanIn
  [n,m] = size(Latmat);
  Ue    = Ue(:,1:length(lon_coads)-1)+diff(Ue,1,2)/2;
  Ue    = [zeros(length(lat_coads),1),Ue,zeros(length(lat_coads),1)];
  q     = find(isnan(Ue));
  Ue(q) = 0;
  
  Ve = Ve(1:length(lat_coads)-1,:)+diff(Ve,1,1)/2;
  Ve = [zeros(1,length(lon_coads));Ve; ...
	    zeros(1,length(lon_coads))];
  q  = find(isnan(Ve));Ve(q)=0;

  % create interpolated fields of salam and temp
  % on Ue, Ve grid.  fill NaN-ed values, or
  % Temp and Salam Ekman convergence/divergence
  % against coasts will be removed.

  Salamx    = [Salam(:,1),Salam,Salam(:,m)];
  Salamx    = fillnans(Salamx,0);
  Salamx    = Salamx(:,1:m+1)+diff(Salamx,1,2)/2;
  Tempx     = [Temp(:,1),Temp,Temp(:,m)];
  Tempx     = fillnans(Tempx,0);
  Tempx     = Tempx(:,1:m+1)+diff(Tempx,1,2)/2;
  landmaskx = [landmask(:,1),landmask,landmask(:,m)];
  landmaskx = landmaskx(:,1:m+1)+diff(landmaskx,1,2)/2;
  landmaskx = ceil(landmaskx);
  Salamy    = [Salam(1,:);Salam;Salam(n,:)];
  Salamy    = fillnans(Salamy,0);
  Salamy    = Salamy(1:n+1,:)+diff(Salamy,1,1)/2;
  Tempy     = [Temp(1,:);Temp;Temp(n,:)];
  Tempy     = fillnans(Tempy,0);
  Tempy     = Tempy(1:n+1,:)+diff(Tempy,1,1)/2;
  landmasky = [landmask(1,:);landmask;landmask(n,:)];
  landmasky = landmasky(1:n+1,:)+diff(landmasky,1,1)/2;
  landmasky = ceil(landmasky);
  
  % replace NaN's with 0 where landmask=1
  % (no Ekman transport over land, which
  % yields correspondingly large conv/divergences)
  q = find(landmaskx==1);Ue(q)=0;
  q = find(landmasky==1);Ve(q)=0;
  clear landmaskx landmasky;

  sdlat = sign(mean(diff(lat_coads)));
  sdlon = sign(mean(diff(lon_coads)));

  % calculate Mek, the Ekman mass convergence (m^3/s)
  Mek(:,:,month) = -sdlon*diff(Ue,1,2)-sdlat*diff(Ve,1,1);
  % calculate Sek, the Ekman salt anomaly conv. ((psu/1e3) m^3/s)
  Sek(:,:,month) = -sdlon*diff(Ue.*Salamx,1,2) -sdlat*diff(Ve.*Salamy,1,1);
  % calculate Tek, the Ekman temperature conv. (deg.C m^3/s)
  Tek(:,:,month) = -sdlon*diff(Ue.*Tempx,1,2) -sdlat*diff(Ve.*Tempy,1,1);

end;    % loop over months

% NaN-out zero values
q = find(Mek==0);Mek(q)=NaN;
q = find(Sek==0);Sek(q)=NaN;
q = find(Tek==0);Tek(q)=NaN;

% divide by number of months, to get average
EkmanInM = EkmanInM/size(sst,1);
EkmanInH = EkmanInH/size(sst,1);
EkmanInS = EkmanInS/size(sst,1);

clear Salamx Salamy Tempx Tempy
clear m n p q rE Taux Tauy Ue Ve Salam Temp Rho
clear Gamma f dlon dlat lat_eq DaysInMonth