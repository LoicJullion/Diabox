function [fbox] = getfluxes(sectfile, properties, reflevel)

%=======================================================================
% GETFLUXES   4.2   93/07/05   Copyright (C) Phil Morgan 1992
%
% [fbox] = getfluxes(sectfile, properties, reflevel)
%
% DESCRIPTION:
%    Calculates the integrated flux of each property for each neutral
%    surface layer at each station pair using the high resolution binned
%    ctd data files.
%
% INPUT:
%    sectfile   = a section filename.  See build_dl.m for details.
%    reflevel   = the "surface" level number to use as the geostrophic
%                 reference level. (old "level of no motion" notion).
%
% OUTPUT:
%    fbox       = integrated flux of each property in layers between 
%                 "surfaces" and "total" flux.    
%
% AUTHOR:   Phil Morgan  20-03-92.
% MODIFIED: Bernadette Sloyan 07-09-95
% MODIFICATION: Now able to specify different reference level for 
%               a single section
% Modification2: Now able to save computed geostrophic velocity using
%                specified reference level (07-10-95)
%=======================================================================
%
% @(#)getfluxes.m  Revision: 4.2  Date: 93/07/05
% @(#)getfluxes.m  Revision: 4.3  No need to enter the names manually
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
%=======================================================================
%% Load variables
% load pressure, number of density surfaces and number of stat pairs
[bin_press, Layer_press]  = loadvars(sectfile,'bin_press','Layer_press');
[deldistkm,nsurfs,npairs] = loadvars(sectfile,'deldistkm','nsurfs','npairs');
nlayers = nsurfs+1; % Number of layers
nbins = length(bin_press); % Number of pressure bins

%% Calculate geostrophic velocity
% Load the baroclinic velocity referenced to the surface and the 
[vel, com_rows, surfPair_bin] = ...
     loadvars(sectfile,'binPair_vel','com_rows','surfPair_bin');

if reflevel == -1
   ref_bin = com_rows;
elseif reflevel == 0
   ref_bin = ones(1,length(com_rows));
else   %   1<=reflevel<=maxlevels
   ref_bin = surfPair_bin(reflevel,:);   
end %if

% Calculate geostrophic velocity
gvel = geovelrel(vel,ref_bin,com_rows,sectfile);

% save gvel to file
load dir_loc.mat gvel_dir
gvelfile = sectfile;
q = find(gvelfile=='/');
gvelfile(1:max(q)) = [];
gvelfile = ['gvel_',gvelfile];
disp(['Saving gvel to file {' gvel_dir '/' gvelfile  '} '])
command = ['save ' gvel_dir '/' gvelfile ' gvel'];
setstr(command);
eval(command);

%% Calculate the fluxes for the individual properties
[nproperties,~] = size(properties);
fbox = [];
disp('Calculating Flux Per Depth for Properties...')
for iprop = 1:nproperties
   PropName = [charword( properties(iprop,:) )];
   disp(['      ' PropName])
   if iprop==1
      binPair_prop = ones(nbins,npairs);
   else
      Var          = ['binPair_' PropName];
      Expression   = ['loadvars(sectfile,Var)'];
      command      = ['binPair_prop = ' Expression ';'];
      eval(command)
   end %if
  
   command = ['FperDepth_prop = binPair_prop .* gvel;'];
   eval(command)
   
   for ipair = 1:npairs
      Expression   = ['FperDepth_prop(:,ipair) * deldistkm(ipair) * 1000'];
      command      = ['FperDepth_prop(:,ipair) = ' Expression ';'];
      eval(command)
   end %for

   LayerF_prop   = [];

   for ipair = 1:npairs
     LayerFPair_prop   = NaN*ones(1,nlayers);
  
     good_p      = Layer_press(:,ipair); 
     good_surf   = find(~isnan(good_p)); 
  
     if ~isempty(good_surf)
         F_good = integrate(bin_press,FperDepth_prop(:,ipair),good_p(good_surf)); 
       for igood = 1:length(F_good)
         LayerFPair_prop(good_surf(igood)) = F_good(igood);
       end %for
     end %if
 
     LayerF_prop   = [LayerF_prop LayerFPair_prop'];

   end %for

   LayerF_prop  = change(LayerF_prop, '==',NaN,0);

   fbox = [fbox; LayerF_prop; sum(LayerF_prop)];
end %for
%%
return

%---------------------------------------------------------------------------

