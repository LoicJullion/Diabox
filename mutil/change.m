function new = change(old,relation,flag,value)
%==========================================================================
% CHANGE     Change values in a matrix
% CHANGE   1.3   92/03/25
%
% new = change(old,relation,flag,value)
%
% DESCRIPTION:
%    Changes the 'flag' values in the matrix "old" to the new "value"
%    according to the "relation".
% 
% INPUT:
%    old      = matrix containing values related to "flag"
%               are to be converted to "value"
%    flag     = values related to "flag" then replaced by "value"
%    value    = replacement value
%    relation = string relation e.g. '<', '>', '=='
%
% OUTPUT:
%    new      = matrix "old" with all flagged values changed
%               to "value" (can be returned to same matrix "old")
%
% EXAMPLE:  A = change(A,'==',NaN,0  )
%           B = change(A,'<', -99,Nan) 
%
% CALLER:   general purpose
% CALLEE:   none
%
% AUTHOR:   Phil Morgan 3-02-92

% @(#)change.m   1.3   92/03/25
% 
% Re-created after 2-2-92 hard disk crash.
% Based on flagnan.m - Steve Rintoul  Dec 90
%          alter.m   - Phil  Morgan    Feb 91
%          convert.m - Peter Mcintosh Aug 91
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

% CHECK INPUT ARGUMENTS CALL
if nargin ~= 4
  error('CHANGE.M: Must have 4 input arguments')
end

if (strcmp(relation,'==') | strcmp(relation,'>')) | strcmp(relation,'<')
    % valid relation
  else
    error(['CHANGE.M: Relation {' relation '} not valid'])
end

% BODY
if isnan(flag)
   replace = find(isnan(old));
else
   eval(['replace = find(old',relation,'flag);']);
end
nreplace = length(replace);
new      = old;
if nreplace>0
   new(replace) = value*ones(1,nreplace);
end %if
return