%==========================================================================
% truncatereanalysis.m
%
% truncate the data from the climate reanalysis (ERA interim) fields by 
% dropping the values of lonvec, latvec, sst, etc. corresponding to vector 
% q.  
% All fields are loaded by loaderasyear.m
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
% Recreate a landmask based on the vectorized data. Good test to see if
% your indeces and your data correspond.
landmask=NaN*ones(length(lat_coads),length(lon_coads));
landmask(index)=mean(sst,1);
landmask=isnan(landmask);
% 
% % Test: If the figure does not look like the landmask in your box, then
% % there is an issue with the vectorized data and/or the indeces.
% figure;
% pcolor(double(landmask)); shading flat
% close

% Get data points outside the box
q = find(lonvec>max(lon)+max(abs(diff(lon_coads))) | ...
	lonvec<min(lon)-max(abs(diff(lon_coads))) | ...
	latvec> max(-50, ...
        max(lat)+max(abs(diff(lat_coads)))) | ...
	latvec<min(lat)-max(abs(diff(lat_coads))));

% Get rid of the points outside the box
lonvec(q)    = [];
latvec(q)    = [];
eminusp(:,q) = [];
netheat(:,q) = [];
sst(:,q)     = [];
sss(:,q)     = [];
taux(:,q)    = [];
tauy(:,q)    = [];

% Isolate the lat and long and landmask
q = find(lon_coads<min(lonvec) | lon_coads>max(lonvec));
lon_coads(q)  = [];
landmask(:,q) = [];

q = find(lat_coads<min(latvec) | lat_coads>max(latvec));
landmask(q,:) = [];
lat_coads(q)  = [];

% Recreate the index corresponding to the truncated matrix 
index = NaN*ones(length(latvec),2);
for n=1:length(lat_coads);
  q=find(latvec==lat_coads(n));
  index(q,1)=n;
end;
for n=1:length(lon_coads);
  q=find(lonvec==lon_coads(n));
  index(q,2)=n;
end;
%rindex=index';
index=(index(:,2)-1)*max(index(:,1)) + index(:,1);
index=index';
Latmat=lat_coads*ones(1,length(lon_coads));
Lonmat=ones(length(lat_coads),1)*lon_coads';

