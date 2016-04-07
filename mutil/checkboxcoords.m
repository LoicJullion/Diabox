function checkboxcoords(boxnum)
%=========================================================================
% function checkboxcoords(boxnum);
%
% AUTHORS: Rick Lumpkin, modified by Loic Jullion to make it user friendly.
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

load dir_loc.mat boxcoord_dir
for loopnum=1:length(boxnum);
    
    tboxnum=boxnum(loopnum);
    checkgeo('geo');
    hold on;

    eval(['load ' boxcoord_dir 'boxcoords',num2str(tboxnum)]);
    
    dlon=min(5,.025*range(lon));
    dlat=min(5,.025*range(lat));
    lonmat=min(lon)-1:dlon:max(lon)+1;
    latmat=min(lat)-1:dlat:max(lat)+1;

    lonmat=ones(length(latmat),1)*lonmat;
    latmat=latmat'*ones(1,size(lonmat,2));
    inbox=NaN*ones(size(lonmat));

    

    q=find(isinpoly(lonmat,latmat,lon,lat));
    plot(lonmat(q),latmat(q),'.');
    hold off;

    if (range(lon)>300)
        axis([-180 180 -90 90])
    elseif (range(lon)>90)
        axis([min(min(lon)-10,-100) max(lon)+10 ...
            min(lat)-10 max(lat)+10]);
    else
        axis([min(min(lon))-.1*range(lon) ...
            max(max(lon))+.1*range(lon) ...
            min(min(lat))-.1*range(lat) ...
            max(max(lat))+.1*range(lat)]);
    end;
    set(gca,'dataaspectratio',[1 1 1])
    title(['Box ',num2str(tboxnum)],'fontsize',15);
    pause(.2);
    
end;    % loop over boxnum