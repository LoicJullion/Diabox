function [dist,bearing] = dobox_distance(lat,long,units)

%===================================================================
% DIST   1.3  92/03/23
%
% [dist,bearing] = distance(lat,long, {units} )
%
% DESCRIPTION:
%   Calculate distance between two positions on glode using the "Plane
%   Sailing" method.  Also uses simple geometry to calculate the bearing of
%   the path between position pairs.
% 
% INPUT:
%    lat      = decimal degrees (+ve N, -ve S)
%    long     = decimal degrees (+ve E, -ve W)
%    units    = optional string specifing units of distance
%               'nm'  (default) = nautical miles
%               'km'            = kilometres
%
% OUTPUT:
%    dist     = distance between positions in units
%    bearing  = bearing of path between positions in degrees.
%               ie  N=0, E=90, S=180, W=270
%
% EXAMPLE: 
%   lat  = [-40   -40   -50   -50   -40   -40   -40]
%   long = [160   170   170   160   170   -170  170]
%   [dist,bearing]=distance(lat,long,'nm')
%
%   dis =
%     459.6267  600.0000  385.6726  734.8469  919.2533  919.2533
%   bearing =
%      90.0000  180.0000  270.0000   35.2644   90.0000  270.0000
%
% CALLER:   general purpose
% CALLEE:   none
%
% AUTHOR:   Phil Morgan and Steve Rintoul 10-02-92
%==================================================================

% @(#)distance.m   1.3   92/03/23
% 
% REFERENCE:
%    The PLANE SALING method as descriibed in "CELESTIAL NAVIGATION" 1989 by
%    The Australian Antartic Division (CSIRO lib# B623.89 GOR).
%
%    Around Equator 1 degree (lat)= 60 arc mins = 60 n.miles
%    The n.mile per min longitude decreases from 1 n.mile at equator
%                                           to   0 n.mile at the Pole
%    At intermediate latitudes n.miles per min longitude = cos(lat).
%    This distance is called the "Departure".
%
%    Distances between two positions is calcated by simple geometry using
%    the difference in latitudes in n.miles (DLAT) and the departure (dep)
%    which is calculated from the difference in longitude.
%
%      DLAT = dlat*DEG2NM    = diff in lat in n.miles
%      DLON = dlon*DEG2NM    = diff in lon in n.miles AT equator
%      dep  = DLON * cos(ave_latitude)
%      
%      dist = sqrt( DLAT^2 + dep^2)    in m.miles
%           = sqrt( (dlat*DEG2NM)^2 + (dlon*DEG2NM*cos(ave_lat))^2 )
%           = DEG2NM * sqrt( dlat^2 + (dlon*cos(ave_lat))^2 )
%
% TESTING:
%   dis =
%     459.6267  600.0000  385.6726  734.8469  919.2533  919.2533
%   anglex =
%            0  270.0000  180.0000   54.7356         0  180.0000
%   bearing =
%      90.0000  180.0000  270.0000   35.2644   90.0000  270.0000
%--------------------------------------------------------------------

DEG2RAD = (2*pi/360);
RAD2DEG = 1/DEG2RAD;
DEG2MIN = 60;
DEG2NM  = 60;
NM2KM   = 1.8520;    % Defined in Pond & Pickard p303.

npositions = length(lat);
ind=1:npositions-1;     % index to first of position pairs

dlong = diff(long);
if any(abs(dlong)>180)
   flag = find(abs(dlong)>180);
   for ii=1:length(flag)
     dlong(flag(ii))= -sign(dlong(flag(ii))) * (360 - abs(dlong(flag(ii))) );
   end %for
end %if
latrad = abs(lat*DEG2RAD);
dep    = cos( (latrad(ind+1)+latrad(ind))./2 ) .* dlong;
dlat   = diff(lat);
dist   = DEG2NM*sqrt(dlat.^2 + dep.^2);

if nargin==3 
  if strcmp(units,'km')
    dist = dist * NM2KM;
  end %if
end %if  

% CALCUALTE ANGLE TO X AXIS
phaseangle  = angle(dep+dlat*sqrt(-1))*RAD2DEG;

% CONVERT ANGLE TO X AXIS TO BEARING
% "anglex" GIVES ANGLE TO X AXIS IN DEGREES. 0...+180 = +ve y axis
%                                            0...-180 = -ve y axis
% Correct -ve y axis values s.t.  E=0, N=90, W=180, S=270.
% ie is anticlockwise direction

anglex = phaseangle;

for i=1:length(anglex)
  if anglex(i)<0         % in 3rd/4th quadrants
    anglex(i) = 360+anglex(i);
  end %if
end %for

anglex = anglex - 90;    % set N=0

for i=1:length(anglex)
  if anglex(i)<0         % fix those in 1st quad
    anglex(i)= 360 + anglex(i); 
  elseif anglex(i)==360  % full circle, set 360 to 0
    anglex(i) = 0;
  end %if
end %for

% REQUIRE BEARING s.t. N=0, E=90, S=180, W=270
bearing = 360 - anglex; % To reflect thru N-S s.t clockwise angle rotation
%--------------------------------------------------------------------

