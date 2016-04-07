function [word] = charword(astring)
%===========================================================================
% WORD  4.1   92/11/12   Copyright (C) Phil Morgan 1992
%
% [word] = charword(astring)
%
% DESCRIPTION:
%    Extracts the first word from a string.  A word is delimited by spaces.
%
% INPUT:
%    astring    = a string of characters
%
% OUTPUT:
%    word       = substring of chars that constitute the first word in astring
%
% EXAMPLE:   word = charword('  hello  there')
%            word = 'hello'
%        
%
% CALLER: general purpose
% CALLEE: none
%
% AUTHOR: Phil Morgan  92-10-16
%
% @(#)charword.m  Revision: 4.1   92/11/12   Copyright (C) Phil Morgan 1992
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
ichars =  find( astring~=blanks(1) );
first = ichars(1);
%astring( first:length(astring) );
iendblanks = find( astring( first:length(astring) ) == blanks(1) );
if isempty(iendblanks)
  last = length(astring);
  word = astring(first:length(astring));
else
  word = astring(first:first+iendblanks(1)-2);  
end %if

return
%-------------------------------------------------------------------
