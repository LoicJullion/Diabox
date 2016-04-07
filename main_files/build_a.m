function [A, nlayers, nsectpairs] = build_a(geometry, sectfiles, properties)

%===================================================================
% build_a   4.3  93/07/01   Copyright (C) Phil Morgan 1991
%
% [A, nlayers, nsectpairs] = build_a(geometry, sectfiles, properties)
%
% DESCRIPTION:
%    Constructs the matrix A in the linear system Ax=b from the specified
%    defined by "geometry", "sectfiles" and "properties".
% 
% INPUT:
%    geometry   = DIM( nboxes, nsections )
%                 Describes the geometry(ibox,isect) of how section "isect"
%                 is used to construct box "ibox".
%                 -1 = section "isect" directed clockwise around box "ibox"
%                  0 = section "isect" is NOT used to construct box "ibox"
%                 +1 = section "isect" directed anti-clockwise around box
%
%    sectfiles  = DIM( nsections, * )
%
%    properties =  DIM(nproperties,*)
%
% OUTPUT:
%    A          = matrix A of the linear system Ax=b
%    nlayers    = number of layers in box
%    nsectpairs = array of number of stations pairs in each section
%
% EXAMPLE:  [A, nlayers, nsectpairs] = build_a(geometry, sectfiles, properties)
%
% CALLER:   dobox.m
% CALLEE:   loadvars.m, charword.m 
%
% AUTHOR:   Phil Morgan 91-11-29
%==================================================================
%
% @(#)build_a.m  Revision: 4.3  Date: 93/07/01
%
% @(#)build_a.m  Commented by Loic Jullion - 2015
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

disp ('get started of build_a.m')

% Define number of boxes and number of stations
[nboxes, nsections] = size(geometry);
% Number of properties
[nproperties,m]     = size(properties);

disp(' ')

% Initialise the A matrix that will contain all the terms of the equations
A      = [];
% Create A
for isect = 1:nsections
  sectfile =  charword( sectfiles(isect,:) );  % strip blanks off name
  disp(' ')
  disp(['build_a: section ' num2str(isect) ' from file {' sectfile '}'])


  A_box = [];
  for iprop = 1:nproperties
      propchars = charword( properties(iprop,:) );
      disp(['loading property ' propchars])
      prop = loadvars(sectfile,['Area_' propchars]);
      A_box = [A_box; prop; sum(prop)];
  end %for

  
  [numlayers(isect), nsectpairs(isect) ] = size(prop);
  nlayers = numlayers(1);     % should all be the same!
  if isect > 1
    if numlayers(isect) ~= numlayers(isect-1)
      error(['build_a: number of layers in files' sectfiles(isect,:) ...
         'and' sectfiles(isect-1,:) ' not the same '])
    end
  end
  
  
  A_sect = [];
  for ibox = 1:nboxes
    A_sect = [A_sect; geometry(ibox,isect)*A_box];
  end
  

  A = [A A_sect];
  
end

return
%--------------------------------------------------------------------



