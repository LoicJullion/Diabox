function [Xcols,XCW] = get_diapvel(filename,nboxes,dir)
%==========================================================================
% function  [Xcols,XCW] = get_diapvel(filename,nboxes,dir)
% Load the variables required to estimate interior diapycnal fluxes.
% This is the only compulsory terms (in addition to the geostrophic
% velocity) for the inversion. All other extra unknowns (air-sea fluxes,
% Ekman...) will be added later. Interior diapycnal fluxes are usually key
% to obtain a good inversion and therefore they are added in the main core
% of Diabox.
% 
% For each box, calculate the mean temperature, salinity and area of the
% layer interface (glevels) within the box and save them in a mat file.
%
% INPUT:
%   - geometry: filename: Name of the file you wish to save (or load) the 
%     variables to (from).
%   - nboxes: Number of boxes.
%   - dir: Location of the files where the data for the diapycnal fluxes
%     are stored.
%
% OUTPUT:   
%   - Xcols: extra columns to problem
%   - XCW: extra column weights
%
% The mat files required to calculate the diapycnal fluxes are named:
% mprop_diap_boxnumber.mat.
%
%!!!!!!!!!!!!!!!
% BE CAREFUL THAT THE BOX NUMBER USE IN THE FILE AND THE BOX NUMBER IN THE
% INVERSION AGREE, OTHERWISE YOU WILL USE INFORMATION ABOUT THE INTERIOR OF
% THE WRONG BOX
%!!!!!!!!!!!!!!!
%
% AUTHOR: Loic Jullion (2015) based on Rick Lumpkin xcolint.m 
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

% Check that the directory exist
if exist(dir,'dir') ~= 7
   warning(['Directory ',dir,' does not exist']) 
   return;
end

% Check if the extra column have already been created and saved in a matlab
% file.
filename = [dir filename '.mat'];
if exist(filename,'file') == 2
   % if yes, load the file
   eval(['load ' filename ' XCW Xcols']);
else
%if not, create it and save it
    for i = 1:nboxes 
      command =['mprop_diap',int2str(i), ' area_layer msalt mptemp'];
      
      % load data
      eval(['load ' dir command]);
            
      % define salinity anomaly (!!!This can be changed but it also has to be
      % changed in prepctd.m!!!)
      msaltan= (msalt - 35)*1e-3;

    %... use diag to form the diagonal matrix with areas on the 
    %... diagonal and -area_layer on the next lower diagonal.

      m = length(area_layer);
      [lg,icolg] = lastgood(mptemp','0')

      mcol = diag(area_layer) + diag(-area_layer(1:m-1),-1);

      tcol = diag(area_layer.*mptemp) + diag(-area_layer(1:m-1) .*mptemp(1:m-1), -1);

      scol = diag(area_layer .* msaltan) + diag(-area_layer(1:m-1) .*msaltan(1:m-1), -1);

    %... add two rows of zeros to end of each diagonal matix. These 
    %... balance the xcol to give 18 rows for each property flux.
    %... The last good value of the xcol is the flux across the bottom
    %... boundary. For the Pacific this is at 895 dbars.

      mcol = [mcol; zeros(2,m)];

      tcol = [tcol; zeros(2,m)];  

      scol = [scol; zeros(2,m)];

    %... Form xcol matrix for each box

      xcol0 = zeros(size(mcol));

      eval(['Xcols',int2str(i),' = ', ...
                                    '[mcol xcol0 xcol0;', ...
                                    ' xcol0 scol xcol0;', ...
                                    ' xcol0 xcol0 tcol];']);


    end;		% loop over boxes

    % Form overall Xcols matrix 
    Xcols=Xcols1;
    n=size(Xcols1,1);m=size(Xcols1,2);
    for i=2:nboxes
      eval(['Xcols=[Xcols,zeros((i-1)*n,m);', ...
        'zeros(n,(i-1)*m),Xcols',num2str(i),'];']);
    end;

    Xcol_int = change(Xcols,'==',NaN,0);

    [nr,nc] = size(Xcols);

    XCW_int =  ones(1,nc);

    Xcols = Xcol_int;
    XCW   = XCW_int;

    eval(['save ' filename ' Xcols XCW']); % Save  
end % if
return

