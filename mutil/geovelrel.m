function gvel = geovelrel(vel, ref_rows, last_rows, sectfile)
%======================================================================
% GEOVELREL   1.6  92/11/20   Copyright (C) Phil Morgan 1992
%
% gvel = geovelrel(vel, ref_rows, {last_rows} )
%
% DESCRIPTION:
%    Calculates the geostrophic velocity relative to a "reference" level
%    or the last valid velocity (if ref_row is greater than last_row) given
%    the velocity field .
%
% INPUT:
%    vel       = geostophic velocity relative to sea surface [m/s]
%                matrix(nbins,npairs)
%    ref_rows  = row of "vel" at each pair to set as the 'reference level'.
%    last_rows = optional row number for last valid velocity value to use.
%                        (set to NaN by default)
% 
% OUTPUT:
%    gvel     = geostrophic velocity relative to the reference level.
%               +ve velocities are such that flow is INTO paper when
%                   "light is on the left"  in the S hemisphere OR
%                   "light is on the right" in the N hemisphere               
%
% AUTHOR: Phil Morgan [11-03-92]
% MODIFICATION:
%   Now has option to set non-zero reference velocity
% Modified:
%   Bernadette Sloyan [11-09-95]
%
% @(#)geovelrel.m   Revision: 1.6   92/11/20  (C) Phil Mrgan
% Supercedes "getgvel.m, avegvel.m by Phil Morgan 14-6-91
% AUTHORS: Rick Lumpkin.
%
% MODIFIED: Loic Jullion.
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
if nargin == 2
  last_rows = NaN*ones(ref_rows);
end %if

% ASSIGN THE REFERENCE VELOCITIES TO THAT AT THE ref_depth OR
% IF THE BOTTOM IS BEFORE ref_depth THEN TO BOTTOM VELOCITY.

%disp(' ') 
%disp('------ Reference Velocity Menu -------')
%disp('   1) Set Reference Velocity to be zero       ') 
%disp('   2) Set Non-zero Reference velocity from *.mat file  ')
 
%RV_opt = input([' Select a menu number:  ']);
RV_opt=2;

file_RV=sectfile;
q=find(file_RV=='/');
file_RV(1:max(q))=[];
file_RV=['refvel_',file_RV];
 
if RV_opt == 1 
  for icol = 1:length(ref_rows)
      if isnan(ref_rows(icol))
           refvel(icol)   =  vel( last_rows(icol), icol );
      else
        if isnan( vel(ref_rows(icol),icol) )
           refvel(icol)   =  vel( last_rows(icol), icol );
        else   % ok data
           refvel(icol)   =  vel( ref_rows(icol), icol );
        end %if
      end %if
  end %for
gvel = vel;

elseif RV_opt == 2
  load dir_loc.mat refvel_dir
  eval(['load ',refvel_dir,'/',file_RV]); 
  RV = ones(length(vel),1)*RV;
  gvel = vel + RV;
    for icol = 1:length(ref_rows)
      if isnan(ref_rows(icol))
           refvel(icol)   =  vel( last_rows(icol), icol );
      else
        if isnan( vel(ref_rows(icol),icol) )
           refvel(icol)   =  vel( last_rows(icol), icol );
        else   % ok data
           refvel(icol)   =  vel( ref_rows(icol), icol );
        end %if
      end %if
  end %for
end %if

for icol=1:length(refvel)
   gvel(:,icol) = gvel(:,icol)-refvel(icol);
end %for
return
%==========================================================================
return

