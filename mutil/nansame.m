function [same, mask] = nansame(matrix1, matrix2)
%==========================================================================
% NANSAME   4.1  92/11/12   Copyright (C) Phil Morgan 1992
%
% [same, mask] = nansame(matrix1, matrix2)
%
% DESCRIPTION:
%    Compares the masks of two matrices containing NaNs
%    at the position of bad data and verifies whether the
%    two matrices have NaNs at identical positions. matrix1 - matrix2
% 
% INPUT:
%    matrix1  = matrix to campare NaNs
%    matrix2  = matrix to campare NaNs
%
% OUTPUT:
%    same     = 1 if matrices have NaNs at same positions
%               0 else
%    mask     = matrix mask of NaNs with values of
%                   1   if matrix1 has Nan matrix2 does not
%                   0   if matrix1 == matrix2 NaN positions
%                  -1   if matrix1 has no Nan but matrix2 does
%
% EXAMPLE:  a  = [ 1  2  NaN NaN]
%           b  = [ 1 NaN  3  NaN]
%           [same,mask] = nansame(a,b)
%           same = 0
%           mask = 0    -1     1     0
%
% CALLER:   general purpose
% CALLEE:   none
%
% AUTHOR:   Phil Morgan 92-02-03
%
% @(#)nansame.m  Revision: 4.1    92/11/12   Copyright (C) Phil Morgan 1992  

% comparing masks for matrix1 - matrix2
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
mask1 = isnan(matrix1);
mask2 = isnan(matrix2);
mask  = mask1 - mask2;

badcols = find(any(mask));
n       = length(badcols);
if n>0
  same = 0;
else
  same = 1;
end







