function [out1, out2, out3, out4] = loadvars(file2load, arg1, arg2, arg3, arg4)
%==========================================================================
% LOADVARS   3.2  93/07/05
%
% [out1, out2, out3, out4] = loadvars(file2load, arg1, arg2, arg3, arg4)
%
% DESCRIPTION:
%    Loads named variables from a specified MATLAB file instead of
%    all variables as is the normal default with "load".
% 
% INPUT:
%    file2load = string of MATLAB filename with extension ".mat"
%    arg1     = string of variable name to retrieve
%    arg2     =             "
%    arg3     =             "
%    arg4     =             "
%
% OUTPUT:
%    out1     = corresponding output of input argument
%    out2     =             "
%    out3     =             "
%    out4     =             "
%
% EXAMPLE:  [lat, long] = loadvars('ns89_28S','lat','long')
%
% CALLER:   general purpose
% CALLEE:   none
%
% AUTHOR:   Phil Organ 14-01-92
%
% @(#)loadvars.m   3.2   93/07/05
% 
% Re-created after 2-1-92 hard disk crash.
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
% CHECK CALL AND EXISTENCE OF FILE
if nargin-1 ~= nargout
  error(['LOADVARS.M: Number of variables in & out do not match.'])
else
  % need nout since nargout reset with "load" command later
  nout = nargout;
  nin  = nargin;
end

if ~isstr(file2load)
  setstr(file2load)
end

% if exist(file2load) ~= 2
%   error(['loadvars.m: File {' file2load '} not found'])  
% end

% LOAD ALL VARIABLES AND OUTPUT THOSE REQUIRED
%disp(' ')   % entered as work around for MATLAB bug1.
command = ['load ' file2load];
setstr(command);
eval(command)

if     nout == 1
  out1 = eval(arg1);
elseif nout == 2
  out1 = eval(arg1);
  out2 = eval(arg2);  
elseif nout == 3
  out1 = eval(arg1);
  out2 = eval(arg2);  
  out3 = eval(arg3);  
elseif nout == 4
  out1 = eval(arg1);
  out2 = eval(arg2);  
  out3 = eval(arg3);  
  out4 = eval(arg4);    
else
  error(['LOADVARS.M: Wrong number arguments in = ' nin ' out = ' nout])
end
return