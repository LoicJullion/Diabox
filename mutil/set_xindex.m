function [row_beg,row_end,col_beg,col_end] = ... 
         set_xindex(Xrows,Xcols,row_beg,row_end, ...
                    col_beg,col_end,nboxes,nsections)
% SET_XINDEX   Set index for pointers for extra rows/columns
%==========================================================================
% SET_XINDEX  1.4   93/10/13   Copyright (C)  Phil Morgan  1992
%
% DESCRIPTION:
%    Set index for pointers for extra rows/columns.
%
% AUTHOR:   Phil Morgan 93-07-05
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

if nargin ~= 8
   error('set_xindex.m: wrong number of input parameters')
end %if

if ~isempty(Xrows)
  Xbox = nboxes + 1;
  [mXrows,nXrows]=size(Xrows);
  row_beg( Xbox ) = row_end( nboxes ) + 1;
  row_end( Xbox ) = row_end( nboxes ) + mXrows;  
  clear Xbox mXrows nXrows
end %if

if ~isempty(Xcols)
  Xsect = nsections + 1;
  [mXcols,nXcols]=size(Xcols);
  col_beg( Xsect ) = col_end( nsections ) + 1;
  col_end( Xsect ) = col_end( nsections ) + nXcols;    
  clear mXcols nXcols Xsect
end %if

return
%--------------------------------------------------------------------
