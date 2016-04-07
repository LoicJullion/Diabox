function gvel2 = sw_bottri_gpan(gpan,lon,lat)
%==========================================================================
% function gvel2 = sw_bottri_gpan(gpan,lon,lat)
% 
% Calculate the geostrophic velocities in the bottom triangle by
% extrapolating the geopotential anomaly as in Thompson and Heywood, 2008.
% 
% For each station, geopotential anomaly (gpan) at the shallowest station 
% is extrapolated down to the depth of the deepest station by using the 
% gradient of gpan.
% 
% Input: 
%        - gpan: Geopotential anomaly (calculated using seawater).
%        - lat : Latitudes of the stations.
%        - lon : Longitudes of the stations.
%
% Output: 
%        - gvel2 = Geostrophic velocities extrapolated in the bottom
%        trianges.
%
% --------------
% References:
% Thompson and Heywood, Frontal structure and transport in the north 
% western Weddell Sea, DSRI, Vol. 55, pp. 1229--1251, 2008.
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
% Authors:
%    Loïc Jullion (FSU/MIO).
%
%==========================================================================

% Size.
[press_levels,stations] = size(gpan);

% Initialise counter and variables.
gvel2     = ones(press_levels,stations-1).*NaN;
istatpair = 1;

for istat = 1:stations-1;
    
  % Initialise the new geop anom variable 
  new_phi = ones(press_levels,2).*NaN;

  % Look for the deepest point for 2 adjacent stations
  [phi_dcl1,deepcomm1] = lastgood(gpan(:,istat)); 
  [phi_dcl2,deepcomm2] = lastgood(gpan(:,istat+1)); 
  
  if deepcomm1 < deepcomm2 % If the first station is shallower
      % Calculate the gradient of geop anom.
      phi_dcl = phi_dcl1;
      grad_phi = diff(gpan(deepcomm1:-1:1,istat));
      
      % Filter
      % WindowSize = 4;
      % grad_phi   = filter(ones(1,WindowSize)/WindowSize,1,grad_phi);
      
      % Fill the geop anom down to the deepest point
      new_phi(1:deepcomm1,1) = gpan(1:deepcomm1,istat);
      
      % In shallow stations the bot tri. might be greater than the depth of
      % the station, then do a trick.
      length_bottri = abs(deepcomm2 - deepcomm1);
      
      if length_bottri < deepcomm1 
         for z = deepcomm1+1:deepcomm2
             % Extrapolate down to the deepest point
             new_phi(z,1) = new_phi(z-1,1)-grad_phi(z-deepcomm1); 
         end
      else % If bot tri is greater than total depth, then just interpolate 
           % gradient to create more poins.
         grad_phi = interp1([1:deepcomm1-1],grad_phi,1:0.1:deepcomm1-1); 
         for z = deepcomm1+1:deepcomm2
             % Extrapolate down to the deepest point
             new_phi(z,1) = new_phi(z-1,1)-grad_phi(z-deepcomm1);
         end 
      end
      % Now fill in the other station
      new_phi(:,2) = gpan(:,istat+1);

  elseif deepcomm2 < deepcomm1 % If the second station is shallower
      % Calculate the gradient
      phi_dcl = phi_dcl2;
      grad_phi = diff(gpan(deepcomm2:-1:1,istat+1));
      new_phi(1:deepcomm2,2) = gpan(1:deepcomm2,istat+1);
      
      length_bottri = abs(deepcomm2 - deepcomm1);
      
      if length_bottri < deepcomm2
         for z = deepcomm2+1:deepcomm1
             new_phi(z,2) = new_phi(z-1,2)-grad_phi(z-deepcomm2);
         end
      else
         grad_phi = interp1([1:deepcomm2-1],grad_phi,1:0.1:deepcomm2-1); 
         for z = deepcomm2+1:deepcomm1
             new_phi(z,2) = new_phi(z-1,2)-grad_phi(z-deepcomm2);
         end 
      end
      new_phi(:,1) = gpan(:,istat);
  elseif deepcomm1 == deepcomm2 % If both have the same depth
      new_phi(:,1:2) = gpan(:,istat:istat+1); % No bottom triangle.
  end 
  
  % Now calculate the geostrophic velocity for this station pair
  gvel2(:,istatpair) = sw_gvel(new_phi(:,1:2),lat(istat:istat+1),lon(istat:istat+1));
  
  % Move to the next station pair.
  istatpair = istatpair+1;
  
  clear grad_phi deepcomm1 deepcomm2 phi_dcl2 phi_dcl1 new_phi
end