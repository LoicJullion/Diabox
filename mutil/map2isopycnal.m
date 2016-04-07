function [prop_iso] = map2isopycnal(prop,dens_vec,isopycnals)
%==========================================================================
% function [prop_iso] = map2isopycnal(prop,dens_axis,isopycnals)
%
% Map the profile of a given property (prop) onto a given set of
% isopycnals. This is a relatively simple way of mapping onto isopycnals in
% 1D.
%
% INPUT:
%      - prop: Property (T, S, O2...) profile to be mapped onto isopycnals.
%      - dens_axis: Density profile corresponding to the property profile.
%      - isopycnals: Isopycnals onto which the property will be mappedl
%
% AUTHORS: Loic Jullion.
%
% DATE: July 2015
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

% Pre-allocate memory
prop_iso = ones(length(isopycnals),1).*NaN;

% ensure climatological gamma profiles are monotonic
for ii=2:length(find(isfinite(dens_vec)))
    if dens_vec(ii)<=dens_vec(ii-1)
       dens_vec(ii)=dens_vec(ii-1)+0.00005;
    end % if
end % for ii
        
ind = find(isnan(dens_vec)==0); % Find values that are not NaNs
if length(ind) >2
   prop_iso = interp1(dens_vec(ind),prop(ind),isopycnals);
end
clear ind dens_vec

return