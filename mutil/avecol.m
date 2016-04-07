%------------------------------------------------------
% AVECOL  1.1  92/04/02  Copyright (C) Steve Rintoul 1990
%
% [abar]=avecol(a)
%
% DESCRIPTION:
%    Averages the columns(pairwise) of a matrix
%
% INPUT:
%   a    = matrix(m*n)
%
% OUTPUT:
%   abar = matrix(m*n-1) columns averaged across rows
%
%  @(#)avecol.m   1.1   92/04/02
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
function [abar]=avecol(a)

[m,n]=size(a);
abar=(a(:,2:n)+a(:,1:n-1))/2;
return