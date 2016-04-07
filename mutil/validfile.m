function [fullname,FileOK] = validfile(filename,extension)
%==========================================================================
% VALIDFILE  4.2   92/12/14   Copyright (C) Phil Morgan 1992
%
% [fullname,FileOK] = validfile(filename, {extension} )
%
% DESCRIPTION:
%    Verifies that a file exists else returns flag FileOK = 0. 
%    User can test FileOK to perform  appropriate action.
%
% INPUT:
%    filename  = string that contains a name of a file
%    extension = optional string that contains the extension to a file name
%
% OUTPUT:
%    fullname  = full name of file = [filename extension]
%    FileOK    = 0 if file does not exit, >0 if file exists
%
% EXAMPLE:  [fname,fileok] = validfile(filename,'.mat')
%        
% CALLER: general purpose
% CALLEE: none
%
% AUTHOR: Phil Morgan  92-10-16
%
% @(#)validfile.m  4.2   92/12/14   Copyright (C) Phil Morgan 1992
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
if nargin == 1
  extension = [];
end %if

  % CHECK IF FILE "filename" EXISTS.
  nchars = length(filename);
  nexten = length(extension);

  if nchars > nexten  % may already have extension in filename 
    if ~strcmp(filename(nchars-nexten+1:nchars),extension)
        fullname = [filename extension];
    else
        fullname = filename;      
    end %if
  else   % filename is shorter than extension.  Thus add extension
	fullname = [filename extension];
  end %if

  FileOK = exist(fullname);
 
return

