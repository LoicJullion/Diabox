%==========================================================================
% makeairsea.m
%
% Calculate the surface forcing for transformation of mass, heat and salt 
% in a given region of the ocean.
% In the Lumpkin and Speer inverse model, there was a choice of 5
% climatologies. Here I have updated the code to use the ERA-Interim
% climatologies which seem to be the best in the Southern Ocean.
%
% Variables saved:
%
%       - glayers: mean gamma of the model layers.
%       - glevels: gamma interfaces between the layers.
%       - gamma2: glevels surrounding by 2 boundary gamma (calculated in
%       glev.m).
%	    - Mm: Water mass formation rates from A/S forcing.
%       - Fh: Net heat flux in the model layers.
%       - Ffw: Net freshwater flux in the model layers including river
%       runoff.
%       - Ffw_norivers: Same as Ffw but without the river runoff.
%       - Mmek: The net mass transport across density classes caused by 
%       direct wind forcing.
%       - Mhek: The net heat transport across density classes caused by 
%       direct wind forcing.
%       - Msek: The net salt transport across density classes caused by 
%       direct wind forcing.
%       - flag_fw: Set to zero for no net E-P-R.
%       - flag_h: Set to zero for no net heat.
%       - R: River runoff.
%       - EkmanInM: Net Ekman mass transport into the box as a function of 
%       gamma (in bins defined by glevels).
%       - EkmanInH: Net Ekman heat anom transport into the box as a  
%       function of gamma (in bins defined by glevels). 
%       - EkmanInS: Net Ekman sal anom transport into the box as a function of 
%       of gamma (in bins defined by glevels). 
%       - EkmanInSalt: Net Ekman driven salt flux.
%       - SSSbar: Mean surface salinity of outcropping layer.
%       - SSTbar: Mean surface temperature of outcropping layer. 
%       - SSArea: Mean surface area of outcropping layer.
%       - Fm: Density flux.
%       - si_flux: Silicate flux.
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
clear all;clc;

% input the years from the UWM/COADS flux fields
% to apply to the inverse model.
% e.g. year=[88 89] for Jan. 1988 - Dec. 1989
% year=[-99] for the climatology
year=input('Year(s) (-99 for monthly climatology)? ');
% q=find(year>1900);
% if (length(q)) year(q)=year(q)-1900; end;

flag_fw        = 1;	 % set to zero for no net E-P-R
flag_h         = 1;	 % set to zero for no net heat
error_multfact = .2; % fractional error in air-sea fluxes
load dir_loc.mat riverdir climdir clim_asdir as_dir boxcoord_dir
disp('Loading geometry file "geo".');
file_geom='geo';
eval(file_geom);
nboxes=size(geometry,1);

error_corrdist=20;	% length scale (in grid points)
			% of errors in air-sea fields

%... Load neutral density surface of inverse model
glev;

%%%%%%%%%%
% LOOP THROUGH BOXES
%%%%%%%%%%

for boxnumber=1:nboxes;
clf;
fprintf('\n     BOX %d\n\n',boxnumber);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD AIR-SEA FLUX DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%

 % ECMWF grid: 0.5 by 0.5 deg
sss=NaN*ones(12*length(year),172161);
sst=sss;netheat=sss;eminusp=sss;precip=sss;evap=sss;
taux=sss;tauy=sss;

for y=1:length(year);

    [lon_coads,lat_coads,lonvec,latvec, sss(12*(y-1)+1:12*y,:), index, ...
     sst(12*(y-1)+1:12*y,:), netheat(12*(y-1)+1:12*y,:), ...
     latent(12*(y-1)+1:12*y,:),sensible(12*(y-1)+1:12*y,:),...
     longwave(12*(y-1)+1:12*y,:),shortwave(12*(y-1)+1:12*y,:), ...
     eminusp(12*(y-1)+1:12*y,:),evap(12*(y-1)+1:12*y,:),...
     precip(12*(y-1)+1:12*y,:), taux(12*(y-1)+1:12*y,:), ...
     tauy(12*(y-1)+1:12*y,:)]= loaderayear(year(y)); 
 
end;

%% Remove data outside box using the polygons created for each box (makeboxcoord.m).
eval(['load ' boxcoord_dir 'boxcoords',num2str(boxnumber)]);
truncatereanal;
q=~isinpoly(lonvec,latvec,lon,lat);
netheat(:,q)    = NaN;
eminusp(:,q)    = NaN;
sss(:,q)        = NaN;
sst(:,q)        = NaN;
taux_inbox      = taux;
tauy_inbox      = tauy;
taux_inbox(:,q) = NaN;
tauy_inbox(:,q) = NaN;
%% Load river run off data.
eval(['load ',riverdir,'river_runoff.mat;']);

clear q;
q = find(isinpoly(lon_river,lat_river,lon,lat)==0);
q = [q;find(isnan(lon_river))];
% For now
lon_river(q) = [];
lat_river(q) = [];
river(q,:)   = [];
R(q)         = [];
si_flux(q)   = [];

% calculate net silica input from rivers
si_flux=sum(si_flux);

% assign mean river fw flux to closest eminusp
orig_eminusp=eminusp;
disp('Adding river runoff.');
for riv=1:length(lon_river);
    % Calculate distance between river mouth and all the points
    clear del
    del = sw_dist([latvec lat_river(riv)],[lonvec lon_river(riv)]);
    % Find the closest
    del=find(del==min(del));
    % if more than one point is at the same distance, take the first one.
    % it won't make a big difference since it is then integrated over the
    % whole box. 
    if length(del)>1
        del = del(1);
    end
    % area of the box
    area           = sw_area(lon_coads,lat_coads);
    area           = area(index);
    area           = area(del);
    eminusp(:,del) = eminusp(:,del)-R(riv)/area;
    % R (m^3/s) divided by area (m^2) gives the
    % freshwater flux in units of m/s, same as eminusp.
end;
clear riv del area
% 
if ~flag_fw eminusp = zeros(size(eminusp)); end;
if ~flag_h netheat  = zeros(size(netheat)); end;

%--------------------------------------------------------
% %% Calculate Ekman depth 
% %  (use a mean value of 1.3e-3 for the drag coefficient
% %  in Tau=Rho_air * Cd * U_10^2; Rossby and Montgomery
% %  (1935) value for Ekman depth z_ek=2.38 u_10/sin(lat).
% z_ek = NaN*ones(length(lat_coads),length(lon_coads),12);
% for mon=1:12
%   taux2         = NaN*ones(length(lat_coads),length(lon_coads));
%   tauy2         = taux2;
%   taux2(index)  = taux(mon,:);
%   tauy2(index)  = tauy(mon,:);
%   U10           = sqrt(  sqrt(taux2.^2 + tauy2.^2) / 1.25 / 1.3e-3);
%   z_ek(:,:,mon) = 2.38*U10./abs(sin(2*pi*lat_coads*ones(1,length(lon_coads))/360));
%   clear taux2 tauy2 U10;
% end;
% THIS DOES NOT SEEM TO BE USED ANYWHERE ELSE AFTERWARDS.
%--------------------------------------------------------
% clear lon, lat (boxcoordsN vectors)
clear lon lat
%% Plot sea surface salinity field. 
% Good to check that you have the data in the right box
if(0) % Switch to one if you want to use.
    sss2=NaN*ones(length(lat_coads),length(lon_coads));
    sss2(index)=mean(sss,1);
    plotsections;
    freezeColors;
    hold on;
    pcolor(lon_coads,lat_coads,sss2);
    plot(lon_river,lat_river,'*b');
    colormap(jet)
    caxis([min(min(sss2)) max(max(sss2))])
end
% pause(.1);
%%%

% net heat std. error: error_multfact of net heat
% or 30 W/m^2 over error_corrdist^2
%netheaterr=30;
% E minus P error is maximum of:
% 1. error_multfact of net E-P on outcrop area, or
% 2. 1e-8 m/s E-P over error_corrdist^2

%% Calculate density, salt, heat fluxes 
calcflux;

% calculate matrix of area (m^2) for pixels
area = sw_area(lon_coads,lat_coads);

% Calculate the net heat and freshwater fluxes, and the resulting density 
% flux, as a function of gamma, on the glevels grid.
disp('Integrating transformation over time, area.');
Fm              = zeros(2,length(glevels));
Fh              = Fm;
Ffw             = Fm;
Ffw_norivers    = Fm;
FmE             = [];
FhE             = [];
FfwE            = [];
LayE            = [];
Mmek            = zeros(1,length(glayers));
Epsilon         = Mmek; 
numcells        = Mmek;
Mhek            = Mmek;
Msek            = Mmek;
SSArea          = Mmek;
SSTbar          = Mmek; 
SSSbar          = Mmek;
glayers_edge    = glayers(1:length(glayers)-1);
glayers_edge(1) = -999;
glayers_edge    = [glayers_edge;999];
gamma2_edge     = [-999;glevels;999];

DaysInMonth=[31 28.25 31 30 31 30 31 31 30 31 30 31];
for month=1:size(sst,1);
  ssg2                 = NaN*ones(length(lat_coads),length(lon_coads));
  fm2                  = ssg2;
  fmerr2               = ssg2;
  fh2                  = ssg2; 
  fherr2               = ssg2;
  ffw2                 = ssg2;
  ffwerr2              = ssg2;
  ffw_norivers2        = ssg2;
  sss2                 = ssg2;
  sst2                 = ssg2;
  ssg2(index)          = ssg(month,:);
  fm2(index)           = fm(month,:);
  fmerr2(index)        = fmerr(month,:);
  fh2(index)           = fh(month,:);
  fherr2(index)        = fherr(month,:);
  ffw2(index)          = ffw(month,:);
  ffw_norivers2(index) = ffw_norivers(month,:);
  ffwerr2(index)       = ffwerr(month,:);
  sss2(index)          = sss(month,:);
  sst2(index)          = sst(month,:);
  Mek2                 = Mek(:,:,month);
  Sek2                 = Sek(:,:,month);
  Tek2                 = Tek(:,:,month);
  
  scalefact = DaysInMonth(month-12*floor((month-1)/12)) ./ ...
        mean(DaysInMonth);
    
  % loop through for terms defined on glevels grid
  for i=1:length(glevels);
    q=find(ssg2>=glayers_edge(i) & ssg2<glayers_edge(i+1));
    if (length(q))
      Fm(1,i) = Fm(1,i)+scalefact*nansum( fm2(q).*area(q)./ ...
	            (glayers_edge(i+1)-glayers_edge(i)) );
      Fm(2,i) = Fm(2,i)+scalefact*nansum( fmerr2(q).*area(q)./ ...
	            (glayers_edge(i+1)-glayers_edge(i)) );
      FmE     = [FmE;scalefact*fmerr2(q).*area(q)./ ...
                (glayers_edge(i+1)-glayers_edge(i))];
      Fh(1,i)=Fh(1,i)+scalefact*nansum( ...
        fh2(q).*area(q));
      Fh(2,i)=Fh(2,i)+scalefact*nansum( ...
        fherr2(q).*area(q));
      FhE=[FhE;scalefact*fherr2(q).*area(q)];
      Ffw(1,i)=Ffw(1,i)+scalefact*nansum( ...
        ffw2(q).*area(q));
      Ffw(2,i)=Ffw(2,i)+scalefact*nansum( ...
        ffwerr2(q).*area(q));
      Ffw_norivers(1,i)=Ffw_norivers(1,i)+ ...
	scalefact*nansum(ffw_norivers2(q).*area(q));
      Ffw_norivers(2,i)=Ffw(2,i)+scalefact*nansum( ...
        ffwerr2(q).*area(q));
      FfwE=[FfwE;scalefact*ffwerr2(q).*area(q)];
      LayE=[LayE;i*ones(length(q),1)];
    end; % if length(q)
  end; % loop through layers
  
  % loop through for terms defined on glayers grid
  for i=1:length(glayers);
    q=find(ssg2>=gamma2_edge(i) & ...
	 ssg2<gamma2_edge(i+1));
    if (length(q))
      Mmek(i)=Mmek(i)+scalefact*nansum( ...
	[Mek2(q);0]);
      Msek(i)=Msek(i)+scalefact*nansum( ...
        [Sek2(q);0]);
      Mhek(i)=Mhek(i)+scalefact*nansum( ...
        [Tek2(q);0]);
      SSArea(i)=SSArea(i)+scalefact*nansum(area(q));
      numcells(i)=numcells(i)+scalefact*length(q);
      SSSbar(i)=SSSbar(i)+scalefact*nansum(sss2(q));
      SSTbar(i)=SSTbar(i)+scalefact*nansum(sst2(q));
    end;
  end;
end;

clear fm2 fmerr2 fh2 fherr2 ffw2* ffwerr2 
clear ssg2 Mek2 Sek2 Tek2
% All fields have been summed over each
% month.  Get average by dividing by the
% total number of months (size(sst,1)).
Fm           = Fm/size(sst,1);
Fh           = Fh/size(sst,1);
Ffw          = Ffw/size(sst,1);
Ffw_norivers = Ffw_norivers/size(sst,1);
Mmek         = Mmek/size(sst,1);
Msek         = Msek/size(sst,1);
Mhek         = Mhek/size(sst,1);
SSArea       = SSArea/size(sst,1);

% calculate mean SSS and SST by dividing nansum by number of cells/level
q           = find(~numcells);
numcells(q) = NaN;
SSSbar      = SSSbar./numcells;
SSTbar      = SSTbar./numcells;

if (0)
% allow the possible bias on Ffw, Fh, Fm
% to be correlated over length scale
% error_corrdist (grid points)
for i=1:max(LayE);
  q=find(LayE==i);
  if (length(q))
    multfact=error_corrdist/sqrt(length(q)/size(sst,1));
	% length(q) is divided by size(sst,1), # months,
	% to get a mean # of outcropping points.
    if (multfact>1) multfact=1; end;
    Fm(2,i)=multfact*sum(FmE(q)/size(sst,1));
    Fh(2,i)=multfact*sum(FhE(q)/size(sst,1));
    Ffw(2,i)=multfact*sum(FfwE(q)/size(sst,1));
    Ffw_norivers(2,i)=multfact*sum(FfwE(q)/size(sst,1));
  end;
end;
% NOTE: floor of error_multfact re-added to errors below.
end;

% now add zero endpoints to F* terms,
% so they match the gamma2 grid
Fm           = [zeros(2,1),Fm,zeros(2,1)];
Ffw          = [zeros(2,1),Ffw,zeros(2,1)];
Ffw_norivers = [zeros(2,1),Ffw_norivers,zeros(2,1)];
Fh           = [zeros(2,1),Fh,zeros(2,1)];

% smooth the curves, keeping the values of the
% absolute maximum and absolute minimum
smoothit  = 'F2(k)=(.5*F(j,k-1)+F(j,k)+.5*F(j,k+1))/2;';
smoothit2 = 'F2(k)=(.25*F(j,k-2)+.5*F(j,k-1)+F(j,k)';
smoothit2 = [smoothit2,' +.5*F(j,k+1)+.25*F(j,k+2))/2.5;'];
for i=1:4;
  if (i==1) F=Fm;
  elseif (i==2) F = Ffw;
  elseif (i==3) F = Ffw_norivers;
  else F=Fh;
  end;
  q=[min(find(F(1,:)==min(F(1,:)))) ...
	min(find(F(1,:)==max(F(1,:))))];
  for j=1:2;
    F2=F(j,:);
    for k=2:min(q)-1;eval(smoothit);end;
    for k=min(q)+1:max(q)-1;eval(smoothit);end;
    for k=max(q)+2:length(gamma2)-2;eval(smoothit2);end;
    F(j,:)=F2;
  end;
  if (i==1) Fm=F;
  elseif (i==2) Ffw=F;
  elseif (i==3) Ffw_norivers=F;
  else Fh=F;
  end;  
end;

% calculate formation rates from A/S forcing
Mm      = NaN*ones(size(diff(Fm')'));
Mm(1,:) = -diff(Fm(1,:));
Mm(2,:) = sqrt(Fm(2,1:size(Fm,2)-1).^2 + Fm(2,2:size(Fm,2)).^2);
clear mm mma;

if (0)
  % create a floor for the error terms:
  % minimum error is always at least error_multfact
  % the magnitude of the forcing term.
  q=find(Mm(2,:)<error_multfact*abs(Mm(1,:)));
  Mm(2,q)=error_multfact*abs(Mm(1,q));
  q=find(Ffw(2,:)<error_multfact*abs(Ffw(1,:)));
  Ffw(2,q)=error_multfact*abs(Ffw(1,q));
  Ffw_norivers(2,q)=error_multfact*abs(Ffw_norivers(1,q));
  q=find(Fh(2,:)<error_multfact*abs(Fh(1,:)));
  Fh(2,q)=error_multfact*abs(Fh(1,q));
  q=find(Fm(2,:)<error_multfact*abs(Fm(1,:)));
  Fm(2,q)=error_multfact*abs(Fm(1,q));
end;

figure('PaperType','A4','PaperUnits','centimeters',...
    'InvertHardcopy','off','Color',[1 1 1],'PaperPosition',[1 1 28.7 20],...
    'Position',[20 225 1300 705]);
axes('Position',[0.1 0.50 0.8 0.4],'FontSize',[14])
err1_Fm  = [0 (Fm(1,:)+Fm(2,:))/1e6 0];
err2_Fm  = [0 (Fm(1,:)-Fm(2,:))/1e6 0];
err_Fm   = [err1_Fm err2_Fm];
x_err    = [[gamma2(1),gamma2',gamma2(end)] [gamma2(1),gamma2',gamma2(end)]];
fill(x_err, err_Fm, [.9 .9 .9],'LineStyle','none') % In light gray
hold on;
plot(gamma2,Fm(1,:)/1e6,'r');
line([min(gamma2) max(gamma2)],[0 0],'color','k');
set(gca,'xlim',[min(gamma2) max(gamma2)],'XTickLabel',[]);
hold off;
ylabel('Transformation (Sv)','FontSize',[14]);
box on; grid on;
eval(['title(''Box no ',num2str(boxnumber),' '');']);
axes('Position',[0.1 0.08 0.8 0.4],'FontSize',[14])
plot(glayers,Mm(1,:)/1e6,'r',glayers,Mm(1,:)/1e6,'.r');
set(gca,'xlim',[min(gamma2) max(gamma2)]);
xlabel('\gamma','FontSize',[14]);
ylabel('Formation (Sv)','FontSize',[14]);
line([20 29],[0 0],'color','k');
box on; grid on
eval(['print -depsc ' as_dir 'Fm_boxno',num2str(boxnumber),'.eps;']);


% interpolate Fh, Ffw to glayers
% (they are currently on the grid of gamma2, bracketing
% the glayers values).
Mfw(1,:)= interp1(gamma2',Ffw(1,:),glayers','linear');
Mfw(2,:)= interp1(gamma2',Ffw(2,:),glayers','linear');
MH(1,:) = interp1(gamma2',Fh(1,:),glayers','linear');
MH(2,:) = interp1(gamma2',Fh(2,:),glayers','linear');
Ffw     = Mfw;clear Mfw; 
Fh      = MH; clear MH;

Mfw(1,:)     = interp1(gamma2',Ffw_norivers(1,:),glayers','linear');
Mfw(2,:)     = interp1(gamma2',Ffw_norivers(2,:),glayers','linear');
Ffw_norivers = Mfw;clear Mfw;

% calculate net Ekman-driven salt flux (kg/s)
EkmanInSalt=1.027*(35*EkmanInM + 1000*EkmanInS);

filename='airsea';
if ~flag_fw filename=[filename,'_noFW']; end;
if ~flag_h filename=[filename,'_noH']; end;
filename=[filename,num2str(boxnumber)];
disp(['Saving ',filename,'.mat']);
eval(['save ',as_dir,filename,' glayers glevels gamma2', ...
	' Mm Fh Ffw* Mmek Mhek Msek flag* R', ...
	' EkmanInM EkmanInH EkmanInS EkmanInSalt', ...
	' SSSbar SSTbar SSArea Fm si_flux']);

clear netheat netheaterr ssd ssg sss sss2 sst taux tauy;
clear eminusp fh fherr fm fmerr;
clear ffw ffwerr index latvec lonvec landmask Sek Tek;
clear Mek Fh Mhek Mm Mmek Ffw Msek area q;
clear M Mek EkmanIn lat_river lon_river R

end;            % LOOP THROUGH BOXES
% 
% THE NEXT STEP WAS FROM PREVIOUS MODEL WHEN SEVERAL CLIMATOLOGIES WERE
% USED.
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % now assemble the transformation, formation and
% % Ekman terms from the different climatologies
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% clear all;
% glev;
% nclims=4;	% number of climatologies
% disp('Loading geometry file "geo".');
% file_geom='geo';
% eval(file_geom);
% nboxes=size(geometry,1);
% 
% for boxnum=1:nboxes;
% 
%   load airsea_1_1 EkmanInH
%   tEkmanInH=NaN*ones([size(EkmanInH) nclims]);
%   tEkmanInM=tEkmanInH;
%   tEkmanInS=tEkmanInH;
%   tEkmanInSalt=tEkmanInH;
%   tFfw=NaN*ones([nclims length(glayers)]);
%   tFh=tFfw;tFfw_norivers=tFfw;
%   tFm=NaN*ones([nclims length(gamma2)]);
%   tMhek=NaN*ones([nclims length(glayers)]);
%   tMm=tMhek; tMmek=tMhek; tMsek=tMhek;
%   tSSArea=tMhek; tSSSbar=tMhek; tSSTbar=tMhek;
% 
%   for climnum=1:nclims;
%     eval(['load airsea_',num2str(climnum), ...
% 	'_',num2str(boxnum), ...
% 	' EkmanIn* F* M* SS* si_flux']);
%     tEkmanInH(:,:,climnum)=EkmanInH;
%     tEkmanInM(:,:,climnum)=EkmanInM;
%     tEkmanInS(:,:,climnum)=EkmanInS;
%     tEkmanInSalt(:,:,climnum)=EkmanInSalt;
%     tFfw(climnum,:)=Ffw(1,:);
%     tFfw_norivers(climnum,:)=Ffw_norivers(1,:);
%     tFh(climnum,:)=Fh(1,:);
%     tFm(climnum,:)=Fm(1,:);
%     tMhek(climnum,:)=Mhek(1,:);
%     tMm(climnum,:)=Mm(1,:);
%     tMmek(climnum,:)=Mmek(1,:);
%     tMsek(climnum,:)=Msek(1,:);
%     tSSArea(climnum,:)=SSArea;
%     tSSSbar(climnum,:)=SSSbar;
%     tSSTbar(climnum,:)=SSTbar;
%   end;
%   EkmanInH=mean(tEkmanInH,3);
%   EkmanInM=mean(tEkmanInM,3);
%   EkmanInS=mean(tEkmanInS,3);
%   EkmanInSalt=mean(tEkmanInSalt,3);
%   EkmanInM_range=sum(tEkmanInM,1);
%   EkmanInM_range=reshape(EkmanInM_range, ...
% 	[size(tEkmanInM,2) nclims]);
%   EkmanInM_range=range(EkmanInM_range')';
%   Ffw=[mean(tFfw);range(tFfw)/2];
%   Ffw_norivers=[mean(tFfw_norivers); ...
% 	range(tFfw_norivers)/2];
%   Fh=[mean(tFh);range(tFh)/2];
%   Fm=[mean(tFm);range(tFm)/2];
%   Mhek=mean(tMhek);
%   Mm=[mean(tMm);range(tMm)/2];
%   Mmek=mean(tMmek);
%   Msek=mean(tMsek);
%   SSArea=nanmean(tSSArea);
%   SSSbar=nanmean(tSSSbar);
%   SSTbar=nanmean(tSSTbar);
%   eval(['save airsea',num2str(boxnum), ...
% 	' Ekman* F* M* tF* SS* si_flux' ...
% 	' glevels glayers gamma2']);
% end;	% loop through boxes
% 

