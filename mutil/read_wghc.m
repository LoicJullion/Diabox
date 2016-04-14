function [varargout] = read_wghc(varargin)
%==========================================================================
% Function to read the WOCE Global Hydrographic Climatology (Gouretski and
% Kolterman). 
% 
% INPUTS:
%       - nprop: variable loaded in geo.m. It contains the different 
%         properties used by the inversion. 
%         By default, properties=['mass '; 'salam';'potmp']; 
%
% The netcdf file 'wghc_params.nc' contains 16 variables. By netcdf ID
% order, they are:
% VarID   VarName
%   0      LON
%   1      LAT
%   2      ZAX
%   3      BOT_DEP
%   4      PRES
%   5      TEMP
%   6      TPOTEN
%   7      SALINITY
%   8      OXYGEN
%   9      SILICATE
%  10      NITRATE
%  11      PHOSPH
%  12      GAMMAN
%  13      SIG0
%  14      SIG2
%  15      SIG4
%
%
%==========================================================================

disp('Reading wghc_params.nc');
load dir_loc.mat climdir 
eval([' ncid = netcdf.open(''' climdir 'wghc_params.nc''); ']); % Open the climatology (netcdf)

% load variables.
for i = 1:nargin
    eval(['varID = netcdf.inqVarID(ncid,''' varargin{i} '''); ']);
    temp  = netcdf.getVar(ncid,varID,'double');
    
    if varID ==0 % First variable is longitude. Need to remap onto -180 180 grid
        temp = lon2lon180(temp);
        % Now rearrange the data.
        indneg = find(temp<0);
        indpos = find(temp>=0);

        temp = [temp(indneg);temp(indpos)];
    end
    if varID > 5
        temp = cat(1,temp(indneg,:,:),temp(indpos,:,:));
        temp(temp ==-9) = NaN; % replace bad data with NaNs;
    end
    varargout{i} = temp;
end
netcdf.close(ncid);

disp('Done');

return