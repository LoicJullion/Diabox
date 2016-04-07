function [row_beg, row_end, col_beg, col_end, eqn_beg, eqn_end] = ...
          set_index(nboxes,nsections,nproperties,nsectpairs,nlayers,nboxeqns)

% SET_INDEX   Set various indeces after Abase built.
%==========================================================================
% SET_INDEX  4.5   93/10/12   Copyright (C)  Phil Morgan  1992
%
% DESCRIPTION:
%    Set index for pointers.
% 
% CALLER:   dobox.m
% CALLEE:   none
%
%
% AUTHOR:   Phil Morgan 92-08-03
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
if nargin ~= 6
  error('set_index.m: wrong number of input arguments')
end %if
 
row_end = [];
row_beg = [];
for ibox=1:nboxes
  row_beg = [row_beg ((ibox-1)*nboxeqns + 1) ];
  row_end = [row_end ibox*nboxeqns ];
end
col_beg = [];
col_end = [];
col_beg = cumsum([1 nsectpairs( 1:length(nsectpairs)-1 ) ]);
col_end = cumsum(nsectpairs);

eqn_beg = [];
eqn_end = [];
for iprop = 1:nproperties
   eqn_beg(iprop) = (iprop-1)*(nlayers+1);  % +1 for total eqns
   eqn_end(iprop) = eqn_beg(iprop)+nlayers-1;
end %for

return

%--------------------------------------------------------------------
