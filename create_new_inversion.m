function create_new_inversion(dirname)
%==========================================================================
% Set up a new inversion using DIABOX.
% DIABOX is a modification of DOBOX 4.2, a Matlab software package for 
% constructing, inverting and analyzing box inverse models developed at the 
% CSIRO Climate Change Re- search Program (Morgan, 1994). 
% DIABOX development began in 2000, and has continued under NSF funding 
% (PIs R. Lumpkin and K. Speer, Florida State University). 
% Current development is handled by Rick Lumpkin 
% (Cooperative Institute for Marine and Atmospheric Sciences, Univ. Miami) 
% and Loïc Jullion (FSU, now at MIO in France). 
% DIABOX adds explicit air-sea transformation for outcropping layers, 
% includes greater automation than in DOBOX for often-repeated steps, 
% collects global choices for hydrographic sections into a single file, 
% solves for unknowns via Gauss-Markov inversion, and fixes several bugs 
% identified in DOBOX 4.2 code.
%-----------------
% How it works:
%     1) cd in the directory where you want to create the directory  
%     2) run create_new_inversion: create_new_inversion('test_inversion');
%     3) Copy templates files
% This script will:
%     1) Check if the directory exists and if not, it creates the directory.
%     2) Create the subdirectories:
%           - sections: Hydrographic sections are stored here.
%           - airseafluxes: Results from airseafluxes.m
%           - gvel: Geostrophic velocities.
%           - refvel: A priori reference velocities
%           - output: Results from the inverse model
%           - diapvel: Area, potemp and sal anom of the layer interfaces
%                      inside each box.
%           - boxcoord: Coordinates of a polygon around the rim each box.
%           - xconstraints: Extra constraints that can be applied to the
%                           model
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
%    Loïc Jullion (FSU/MIO). Strongly inspired by Kevin Speer, Rick Lumpkin
%    and the many collaborators who helped me improve my understanding of
%    inverse modelling.
%
%==========================================================================
% Find the path to the directory
s = what('Diabox_v2.0');
diabox_path = s.path;
clear s;
%
if exist(dirname,'dir') == 7
   warning(['Directory ',dirname,' already exists. Delete the existing directory or chose another name']) 
   return;
else
   dirpath = pwd;
   disp(['Creating the directory ',dirname,' in ',dirpath]);
   eval(['! mkdir ',dirname,';']);
end
pause(0.5);
disp(['Creating subdirectorie in ',dirpath]);
eval(['! mkdir ',dirname,'/sections;']);
eval(['! mkdir ',dirname,'/airseafluxes;']);
eval(['! mkdir ',dirname,'/gvel;']);
eval(['! mkdir ',dirname,'/output;']);
eval(['! mkdir ',dirname,'/refvel;']);
eval(['! mkdir ',dirname,'/diapvel;']);
eval(['! mkdir ',dirname,'/boxcoord;']);
eval(['! mkdir ',dirname,'/xconstraints;']);

sect_dir     = ['sections/']; % Directory where sections are saved
as_dir       = ['airseafluxes/']; % Air sea fluxes
gvel_dir     = ['gvel/'];         % Geostrophic vel
refvel_dir   = ['refvel/']; % Reference vel
output_dir   = ['output/']; % Output from the inversion
diapvel_dir  = ['diapvel/']; % Output from the inversion
boxcoord_dir = ['boxcoord/']; % Coordinates of polygons along box's rim
xconst_dir   = ['xconstraints/']; % Extra constraints

riverdir     = '~/Work/Diabox_v2.0_beta/Anc_data/Riverflow/'; % Path to the river flow data
climdir      = '~/Work/Diabox_v2.0_beta/Anc_data/Ocean_clim/'; % Path to ocean reanalysis used to calculate SSS
clim_asdir   = '~/Work/Diabox_v2.0_beta/Anc_data/ERA/'; % Path to climate reanalysis used to calculate air-sea fluxes
bathy_dir    = '~/Work/Diabox_v2.0_beta/Anc_data/bathy/';

% Copy files which will contain inversion specific parameters.
% Geometry file
eval(['! cp ' diabox_path '/geo_template.m ' dirname '/geo.m']);
% Reference velocity file
eval(['! cp ' diabox_path '/uniquerefvels_template.m ' dirname '/uniquerefvels.m']);
% Interior diapycnal velocity file
eval(['! cp ' diabox_path '/xcolint_template.m ' dirname '/xcolint.m']);
% Interior diapycnal velocity file
eval(['! cp ' diabox_path '/add_xconstraints_template.m ' dirname '/add_xconstraints.m']);

eval(['save ',dirname,'/dir_loc.mat *_dir riverdir climdir clim_asdir bathy_dir']);
return;