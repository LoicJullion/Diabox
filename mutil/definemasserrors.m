function [masssecterror,masstotalerror]=definemasserrors(sectfiles,geometry)
%==========================================================================
% function [masssecterror,masstotalerror]= definemasserrors(sectfiles,geometry);
%
% code to calculate a mass error for a section (m/s) due primarily to 
% hydrographic variabilty.
% Factors: close to equator = large, far=small. long=large, short=small.
%
% AUTHORS: Rick Lumpkin.
%
% MODIFIED: Loic Jullion.
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
distl=NaN*ones(1,size(sectfiles,1));
mlatt=distl;
for ix=1:size(sectfiles,1);
    raw_filename = charword(sectfiles(ix,:));
    eval(['load ',raw_filename,'_raw lon lat']);
    distl(ix)=sum(sw_dist(lat,lon,'km'));
    mlatt(ix)=abs(mean(lat));
end;

% first error: large near equator
er1=cot(mlatt*2*pi/360)/2;

% second error: proportional to distance
er2=(distl/5e3);
    
masssecterror=1e6*sqrt(er1.^2+er2.^2);
clear ix distl mlatt er1 er2

% net volume imbalance per box, assuming non-synoptic
% error from section to section are independent.
masstotalerror=sqrt(sum((abs(geometry).* ...
        (ones(size(geometry,1),1)*masssecterror.^2))'));
    return