function [DL,D,bbase,reflevel] = ...
 rhs_menu(geometry,sectfiles,properties,DL,D,bbase,reflevel)
% RHS_MENU  Menu for constructing RHS matrices D, bbase
%===================================================================
% RHS_MENU  4.7   93/10/12   Copyright (C)  Phil Morgan  1992
%
% DESCRIPTION:
%    Script for menu of the constructing RHS
% 
% CALLER:   dobox.m
% CALLEE:   build_dl.m build_e.m
%
%
% AUTHOR:   Phil Morgan 92-08-03
% Modified: Loic Jullion - 2015: Now only include the D variable. I got rid
% of Ekman that is now calculated in the airsea fluxes (airseafluxes.m).
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
%==================================================================

if     nargin == 4
  % ok - minimum allowable input parameters.
  DL = []; D = []; bbase = [];
elseif nargin == 7
  % ok - pass thru any previous values for output parameters DL,E,D,bbase
else
  error('rhs_menu.m: wrong number of input parameters')
end %if

rhs_opt_stop = 0;

rhs_opt = 1;


if isempty(DL)
   disp(' ')
   disp('WARNING: RHS not constructed, reflevel not set.')
   disp(' ')
   disp('Hit a key to continue and construct RHS')
   pause
else
   D = -DL;
   bbase = sum(D')';
end %if

%%%% reflevel now set in geo.m %%%%
[DL] = build_dl(geometry, sectfiles, properties, reflevel);
   
D = -DL;
bbase = sum(D')';
%%%%%

return
%--------------------------------------------------------------------
