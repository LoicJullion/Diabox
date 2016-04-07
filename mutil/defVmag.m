function Vmag=defVmag(section,reflevel,floorVmag,WBCwidth,WBCmag)
%==========================================================================
% function Vmag=defVmag(section,reflevel,floorVmag,WBCwidth,WBCmag);
%
% Define Vmag, the magnitude of the a priori std dev.
% of Vref.  This is the adjustment to the thermal wind
% reference level calculated in the inversion.  
% Assume that the adjustment can be +-250 dbar
% (bottom to bottom+250, if reflevel shoals;
% surface to 250dbar, if reflevel outcrops)
% Use "floorVmag" to define a minimum value for Vmag.
% REFERENCE: Ganachaud (dissertation, 1998, pg. 61).
%
% If WBCmag>0, find all station pairs within WBCwidth
% (km) of the western boundary and set Vmag=WBCmag for
% those pairs.
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
maxvel=.15;

errorrange=250;	% error is +-errorrange (dbar)

eval(['load ',section, ...
	' binPair_vel bin_press surf_press lon lat bottom']);

%f=2*7.2921e-05*sin(lat*2*pi/360);
% calculate Vbt=L^3 * grad(Hbar)
%L=10e3;
%dx=sw_dist(lat,lon,'km')*1e3;
% calculate H (smoothed bottom)
%H=gausssmooth(bottom,[0 cumsum(dx)]',50e3)';
%f=sw_f(lat);
%Vbt=L.^3.*abs(diff(f./H))./dx;
%q=find(Vbt>maxvel);Vbt(q)=maxvel;
%clear L dx H f

% convert lon, lat (stations) to lon, lat (station pairs)
lon=(lon(1:length(lon)-1)+diff(lon)/2)';
lat=(lat(1:length(lat)-1)+diff(lat)/2)';

Vmag=NaN*ones(1,size(binPair_vel,2));
disp('Not doing the topo-varying vref scaling')
Vbt=.01*ones(size(Vmag));

% calculate Vmag for +-errorrange adjustment
for stnp=1:length(Vmag);
  q=min(find(bin_press>=nanmean( ...
	surf_press(reflevel,stnp:stnp+1))));
  if (length(q));
    q=find(bin_press>=bin_press(q)-errorrange & ...
	bin_press<=bin_press(q)+errorrange);
  elseif (isfinite( max(surf_press(1:reflevel,stnp))  ))
	% ref level is "below" bottom, i.e. 
	% bottom ref level for this station pair
    q=max(find(isfinite(binPair_vel(:,stnp))));
    q=find(bin_press>=bin_press(q)-errorrange);
  else	% ref level outcrops, i.e. surf_press
	% is only defined for lighter gamma.
	%    (isfinite( max(surf_press( ...
	%    reflevel+1:nlayers-1,stnp))  )) == 1
    q=find(bin_press<=errorrange);
  end;
  q=binPair_vel(q,stnp);
  q=[min(q) max(q)];
  Vmag(stnp)=abs(diff(q));
end;

% added 3 Oct 2002: put a ceiling on the
% +- errorbar calculation for Vmag
q=find(Vmag>maxvel);
Vmag(q)=maxvel;

q=find(Vmag<Vbt | isnan(Vmag));
Vmag(q)=Vbt(q);
q=find(Vmag<floorVmag);
Vmag(q)=floorVmag;

% impose WBC floor
if (WBCmag>0);
  fprintf('%s : magnifying Vref in WBC.\n',section);
  q=find(lon==min(lon));
  del=NaN*ones(size(lon));
  for i=1:length(del);
    del(i)=sw_dist([lat(i) lat(i)],[lon(i) lon(q)],'km');
  end;
  q=find(del<=WBCwidth);
  for wj=1:length(q);
    Vmag(q(wj))=max([Vmag(q(wj)) WBCmag]);
  end;
end;	% if-adding-WBC loop

return
