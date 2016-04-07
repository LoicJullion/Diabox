function [DL] = build_dl(geometry, sectfiles,properties, reflevel)

%===================================================================
% BUILD_DL   4.3   93/07/01   Copyright (C) Phil Morgan 1992
%
% [DL] = build_dl(geometry, sectfiles, properties, reflevel)
%
% DESCRIPTION:
%    Constructs the the reLative integrated property fluxes between 
%    surfaces using the geostrophic velocity relative to a reference level.
%    DL denotes Data flux reLative.
%
%    Ax = b =  sum(DL')  where 
% 
% INPUT:
%    geometry   = DIM( nboxes, nsections )
%    sectfiles  = DIM( nsections, * )
%    properties =  DIM(nproperties,*)
%    reflevel   = "surface" number that serves as the reference level.
%
% OUTPUT:
%    DL         = DIM(DL) = DIM(Abase)
%                 Block matrix of integrated layer property fluxes.
%
% EXAMPLE:  [DL] = build_dl(geometry, sectfiles, properties, reflevel)
%
% CALLER:   dobox.m
% CALLEE:   charword.m, getfluxes.m 
%
% AUTHOR:   Phil Morgan 12-03-92
%==================================================================
%
% @(#)build_dl.m  Revision: 4.3  Date: 93/07/01
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
%--------------------------------------------------------------------

[nboxes, nsections] = size(geometry);
disp(' ')

DL = [];

for isect = 1:nsections
  sectfile = charword( sectfiles(isect,:));
  disp(['build_dl: section ' num2str(isect) ' from file {' sectfile '}'])

  [DL_box] = getfluxes(sectfile, properties, reflevel(isect));

  % FILL IN A FOR ALL BOXES IN THIS SECTION
  DL_sect = [];
  for ibox = 1:nboxes
    % Multiply DL by geometry (+1 or -1)
    DL_sect = [DL_sect; geometry(ibox,isect)*DL_box];
  end
  
  % FILL A FOR EACH SECTION 
  DL = [DL DL_sect];
 
end

return
%--------------------------------------------------------------------



