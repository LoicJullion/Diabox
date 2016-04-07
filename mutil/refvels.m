function refvels(geometry,sectfiles,rv,vmag)
%==========================================================================
% function refvels(geometry,sectfiles,rv,vmag)
%
% Define a reference velocity for all station pairs in a section.  
% This is used by DIABOX during step 3.
% 
% AUTHOR: Rick Lumpkin
%
% modified by Loic Jullion: Transformed in a function and uniquerefvels is
% now moved to the main diabox_menu.m file
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

load dir_loc.mat refvel_dir

% loop through sections;
for i=1:size(geometry,2);

  sectname=sectfiles(i,:);
  Expression   = ['lon = loadvars(sectname,''lon'')'];
  eval(Expression);
  
  RV=rv(i)*ones(length(lon)-1,1)';
  Vmag=vmag(i)*ones(length(lon)-1,1)';
  q=find(sectname=='/');
  sectname(1:last(q))=[];
  q=find(sectname=='.');
  sectname=sectname(1:q-1);
  rvfilename=[refvel_dir,'/refvel_',sectname];
  eval(['save ',rvfilename,' RV Vmag']);

end;
return
