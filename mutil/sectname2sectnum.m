function sectnum=sectname2sectnum(sectname,sectfiles)
%==========================================================================
% function sectnum=sectname2sectnum(sectname,sectfiles);
%
% Given the name of a file (e.g. 'a08'), or files (['a08_west';'a08_cent';'a08_east']);
% identify the corresponding section number for that particular box, as
% defined in geo.m
%
% AUTHORS: Rick Lumpkin
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
if ~exist('sectfiles')
    geo
else
    sf=sectfiles(1,:);
    for i=1:length(sf);
        if (strcmp(sf(i),'/'))
            sect_dir=sf(1:i);
        end
    end
end

sectnum=[];
for i=1:size(sectname,1);   % loop over the sectnames
  isectname=sectname(i,:);
  % strip the ".mat   "
  q=findstr(isectname,'.');
  q=[q,findstr(isectname,' ')];
  if (length(q))
    isectname=isectname(1:min(q)-1);
  end
  for j=1:size(sectfiles,1);    % loop over sectfiles
    jsectname=sectfiles(j,length(sect_dir)+1: ...
        size(sectfiles,2));
      %length(directory)+size(sectname,2));
    q=findstr(jsectname,'.');
    q=[q,findstr(jsectname,' ')];
    if (length(q))
      jsectname=jsectname(1:min(q)-1);
    end

    if (strcmp(jsectname,isectname))
    	sectnum=[sectnum,j];
    end
  end
end    % loop over sectnames