function [sectfile] = dosection(raw_filename, properties)
% DOSECTION   Process raw_sectfile to sectfile for dobox
%===================================================================
% DOSECTION  4.7   93/10/12   Copyright (C)  Phil Morgan 1992.
%
% sectfile = dosection(raw_filename, properties)
%
% test version
%
% DESCRIPTION:
%    Reads a file "raw_filename" containing binned "properties" data and 
%    calculates the average property and property*area in each layer.
%    Processed section saved to new file.
%
% INPUT:
%    raw_filename   = Raw section file name string
%    properties     = List of property names in raw_filename. (string)
%    sectfile       = optional filename for new section file. Else prompted.
%
% OUTPUT:
%    sectfile       = processed section file name
%    The raw variables are saved under the name "sectfile" (section_raw)
%    and the bottom triangles are saved in add_tri (section_tria)
%    
%
% AUTHOR:  Phil Morgan 12-03-92 Modified by Rick Lumpkin, Takamasa 
%          Tsubouchi and Loic Jullion
%
%==========================================================================
%
% @(#)dosection.m   Revision: 4.7   Date: 93/10/12
% @(#)dosection.m   Revision: 4.8   Date: 10/03/2014
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
%--------------------------------------------------------------------

disp('DOSECTION determines where valid layers exist and calculates')
disp('          mean property values in layers at station pairs')
disp(' ')
TRUE  = 1;
FALSE = 0;

% define file names from raw_filename prefix
dir_path     = ['sections/'];
add_tri      = [raw_filename,'_tria'];
sectfile     = [dir_path,raw_filename];

raw_filename          = charword(sectfile);   % strip blanks
[raw_filename,fileok] = validfile(sectfile,'_raw.mat');
if fileok
   disp(['Loading data file {' sectfile '_raw.mat}'])
   disp(' ')
   command=['load ' sectfile,'_raw.mat'];
   clear Ln;
   char(command);
   eval(command);    
else
   error(['File {' sectfile '} not found'])
end %if
clear fileok command

% Define number of isopycnal surfaces, properties, station pairs and layers
[nsurfs, nstations] = size(surf_press);
[nproperties, n ]   = size(properties);
npairs              = nstations - 1;
nlayers             = nsurfs+1;  % not used  - can be derived from nsurfs
clear n

if nproperties < 2
   error('require at least 2 props.  mass and one other')
end %if

disp('Verifying that properties are self-consistent')
PropName2 = charword( properties(2,:) );
bin_prop2 = ['bin_' PropName2];
for iprop = 3:nproperties
   prop_next = ['bin_' charword( properties(iprop,:) )];
   if ~nansame( eval(bin_prop2) , eval(prop_next) );
      error(['     ' bin_prop2 ' and ' prop_next ' not consistent'])
   else
      disp(['     ' bin_prop2 ' and ' prop_next ' ARE consistent'])    
   end %if
end %for
clear prop_next bin_prop2 iprop 

addbottri = 1; % Switch to 0 if you do not want bottom triangles.
for iprop = 2:nproperties
    PropName = [charword( properties(iprop,:) )];
    disp(['      ' PropName])
    bin_prop = ['bin_' PropName];
    VariableName = ['binPair_' PropName];
    command = [VariableName ' = avecol(' bin_prop '); '];
    eval(command)
    triindex   = sum(isfinite(binPair_vel));
    eval(['tt = sum(isfinite(' bin_prop ')); ']);
    tt = round(diff(tt)/2);
    eval(['ntriindex  = sum(isfinite(' VariableName ')); ']);
    [mn,nn] = size(binPair_vel);
    for tri = 1:nn
        if addbottri == 1
           disp('Adding bottom triangle to the geostrophic velocities'); 
           if tt(tri) >= 0
              eval([VariableName '(ntriindex(tri) + 1:triindex(tri),tri) = '...
                    bin_prop '(ntriindex(tri)+1:triindex(tri) ,tri+1); ']);
           elseif tt(tri) < 0
              eval([VariableName '(ntriindex(tri) + 1:triindex(tri),tri) = '...
                      bin_prop '(ntriindex(tri)+1:triindex(tri) ,tri); ']);
           end %(if tt)
        elseif addbottri == 0
              binPair_vel(ntriindex(tri) + 1:triindex(tri),tri) = ...
                 ones(triindex(tri)-ntriindex(tri),1)*NaN;   
        end %if 
    end %for
end %for

% Pressure of the isopycnals at the station pairs
surfPair_press    = avecol( surf_press );

clear command VariableName PropName bin_prop ntriindex triindex mn nn 
disp(' ');
disp('Verifying binPair_vel...')
if ~nansame( eval(['binPair_' PropName2]), binPair_vel )
  disp(['    binPair_vel NOT consistent with ' PropName2])
  keyboard
  %error('STOPPED by dosection.m')
else
  disp(['    binPair_vel is consistent with ' PropName2])
end %if


disp('Calculating deepest common pressures...')
% This is the max pressure of the binPair properties
Press2Row         = [bin_press [1:length(bin_press)]' ];
[values,com_rows] = lastgood( eval(['binPair_' PropName2]),NaN);

pelp1=fliplr(Press2Row);
com_press = interp1(pelp1(:,1),pelp1(:,2),com_rows);
clear values PropName2 pelp1

%%
%=-=-=-=-=-= figutr of SurfPair -=-=-=-=-=-=-=-=-=-=-

%reply = input('Do you want to draw SurfPair figure? (y/n): ','s');
reply = char('y');
if(strcmp(reply,'y'))

mlon=avecol(lon);
figure;
hold on;
plot(mlon,surfPair_press,'r-.'); axis ij
plot(lon,surf_press,'k.');
%plot(mlon,com_press,'k-');
hold off;
title('comparison between surfPair_press (red) and surf_press (black)');
xlabel('longitude(deg)'); ylabel('pressure(dbar)');

end
%%
disp('Finding valid surfaces above deepest common pressures...')
surfPair_press_valid = surfPair_press;
surfPair_bin         = NaN*ones(size(surfPair_press));

for ipair = 1:npairs
  bad = find( surfPair_press(:,ipair) > com_press(ipair) );
  if ~isempty(bad)
     surfPair_press_valid(bad,ipair) = NaN*ones(length(bad),1);
  end %if
  good = find( ~isnan(surfPair_press_valid(:,ipair) ) );
  if ~isempty( good )
     surfPair_bin(good,ipair) = ...
         fix(interp1(Press2Row(:,1),Press2Row(:,2), ...
		surfPair_press_valid(good,ipair)) );
  end %if
end %for
clear good bad ipair

% Create a matrix that contains the top and bottom of each layer of the
% model: The first row will be the surface, then the pressure of the
% different isopycnals defined and finally the bottom depth.
Layer_press=NaN*ones(nsurfs+2,npairs);
Area_layer=NaN*ones(nsurfs+1,npairs);

Layer_press = [NaN*ones(size(surfPair_press_valid(1,:)));
               surfPair_press_valid;
	           NaN*ones(size(surfPair_press_valid(1,:)))];

[val1,itop] = firstgood(surfPair_press_valid,NaN);
[val2,ibot] = lastgood( surfPair_press_valid,NaN);

% (2013.10.28.) modified to deal with region where you do not have any
% surf_press value.
for ipair = 1:npairs  %length(val1)
    % Fill in the pressure of the top of the first layer. We set to the
    % first pressure from bin_press.
    Layer_press( 1,   ipair) = bin_press(1); 
    % Now we need to add a bottom to the last layer.
    if ~isnan(ibot(ipair))
       % 1) There is at least 1 isopycnal present at the station pair. Then 
       % find which one it is (ibot(ipair)), and add the depth of the
       % bottom triangle (com_press(ipair)) as the lower limit of the layer
       Layer_press( ibot(ipair)+2, ipair) = com_press(ipair);
    elseif isnan(ibot(ipair))
       % 2) If no isopycnal is present at this station pair, then add the 
       % depth of the bottom triangle to the bottom of the first layer.
       Layer_press( 2, ipair) = com_press(ipair);
    end %if
end %for

Layer_bin = NaN*ones(size(Layer_press));

for ipair = 1:npairs
   goodi=find( ~isnan(Layer_press(:,ipair)) );
   % length(goodi) = 0 (ok, no surfs))
   if length(goodi) > 1
%      Layer_bin(goodi,ipair) = fix(table1(Press2Row,Layer_press(goodi,ipair)) );
       Layer_bin(goodi,ipair) = fix(Layer_press(goodi,ipair));
   elseif length(goodi) == 1
      error('dosection.m: must have more than 1 surface to define layers')
   end %if
end %for

clear Press2Row goodi

%%
disp('Calculating thickness of Layers ... ')
% elm 29.9.99 - calculate thickness from integrated frac so that
% bottom triangle area is properly incorporated

% Layer_press is integrated, taking into account the bottom triangle to
% calculate the thickness of the layers
thick = NaN*ones(nlayers,npairs);

for ipair = 1:npairs
    good_p      = Layer_press(:,ipair);
    good_surf   = find(~isnan(good_p));
    F_good = [];
    F_good = integrate(bin_press,frac(:,ipair),good_p(good_surf));
    % INSERT INTEGRATE FLUXES INTO VALID LAYERS, OTHERS LEFT AS NaN
    for igood = 1:length(F_good)
        thick(good_surf(igood),ipair)   = F_good(igood);
    end %for
end %for

thick = change(thick,'==',NaN,0);

%reply2 = input('Do you want to draw thickness figure? (y/n): ','s');
reply2 = char('y');
if(strcmp(reply2,'y'))

    nlayer=size(Layer_press,1);
    
    figure;
    % This plot is misleading. Looking at it, one might think that some
    % isopycnals are not intersecting the topography and stopping in the
    % middle of the water column. This is not the case, when an isopycnal
    % seem to disappear (or appear) in the middle of nowhere, it just 
    % means that this isopycnal is not present at the next/preceding
    % station pair and therefore it must be running into the topography
    % between the station pairs. Therefore at the preceding/next station
    % pair, the layer defined by this isopycnal is not present.
    hold on;
    for II=1:nlayer; 
        plot(Layer_press(II,:),'k-','LineWidth',2); 
    end
    
    npair=length(com_press);
    maxpr=max(bin_press)+10;
    
    for JJ=1:npair;
        sline=[JJ,JJ];
        plot(sline,[0,maxpr],'b--');
    end;
        
    plot(com_press,'r--.','LineWidth',1); 

    hold off;
    axis ij;
    title('Layer_press (black) com_press (red) stationpair(blue)');

    figure;
    hold on
    for II=1:nlayers;
        if(II==1);plot(thick(II,:),'r-','LineWidth',2);end;
        if(II==2);plot(thick(II,:),'g-','LineWidth',2);end;
        if(II==3);plot(thick(II,:),'b-','LineWidth',2);end;
        if(II==4);plot(thick(II,:),'m-','LineWidth',2);end;
        if(II==5);plot(thick(II,:),'k-','LineWidth',2);end;
    end
    legend('layer1','layer2','layer3','layer4','layer5',-1);
    grid on;
    hold off;
    axis ij;
    title('Layer thickness for each layer');

end % if
%%
disp('Calculating mean property value in middle of LayerPair for ')

LayerPair_mass = ones(nsurfs+1,npairs);

for iprop = 2:nproperties
   PropName = [charword( properties(iprop,:) )];
   disp(['     ' PropName])
   binPair_property = eval(['binPair_'   PropName]);

   IntLayer_prop  = [];
     
   for ipair = 1:npairs
     IntLayerPair_prop   = NaN*ones(1,nlayers);
  
     %% FIND VALID NEUTRAL SURFACES FOR A STATION PAIR. incl csurf and bottom
      good_p      = Layer_press(:,ipair); 
      good_surf   = find(~isnan(good_p)); 
  
     %% INTEGRATE PROP DOWN A STATION PAIR BETWEEN VALID SURFACES
     if ~isempty(good_surf)
%$$$    y  = binPair_property(:,ipair)
%$$$ 	x  = bin_press
%$$$ 	xi = good_p(good_surf)
	F_good = [];
% elm 29.9.99 scale bin_pair property by frac so that integrated property
% is properly weighted
%        F_good = integrate(bin_press,binPair_property(:,ipair),good_p(good_surf));
        F_good = integrate(bin_press,frac(:,ipair).*binPair_property(:,ipair),good_p(good_surf));
       %% INSERT INTEGRATE FLUXES INTO VALID LAYERS, OTHERS LEFT AS NaN
       for igood = 1:length(F_good)
         IntLayerPair_prop(good_surf(igood))   = F_good(igood);
       end %for
     end %if
     
     
     %% APPEND INTEGRATION OF ONE STATION PAIR FLUXES INTO A SECTION OF FLUXES
     IntLayer_prop   = [IntLayer_prop IntLayerPair_prop'];

   end %for
   clear F_good binPair_property igood good_p good_surf 
   
   % want Intergrated_property between surfaces (ie in layer) divided
   % by layer thickness.  Avoid division by zero error messages.
   bad = find( thick==0);
   goodthick = thick;
   if ~isempty(bad)
      goodthick(bad) = ones(1,length(bad));
   end %if
   Layer_mean     = IntLayer_prop ./ goodthick;
   Layer_mean     = change(Layer_mean,'==',NaN,0);
   if ~isempty(bad)
      Layer_mean(bad) = zeros(1,length(bad));
   end %if
   clear bad goodthick

   VariableName = ['LayerPair_' PropName];
   command = [VariableName ' = Layer_mean;'];
   eval(command)
   
end %for   

% ---- FIGURE -------

dfig = 0;

if (dfig == 1);
    figure;
    pcolor(LayerPair_temp); shading flat; colorbar; axis ij; caxis([-1.5 4]);
    title('LayerPair_temp');

    figure;
    plot(LayerPair_temp,'r.'); axis([0 6 -1.5 7]);
    title('LayerPair_temp');
    
    figure;
    pcolor(LayerPair_sal); shading flat; colorbar; axis ij; caxis([33.8 34.8]);
    title('LayerPair_sal');

    figure;
    plot(LayerPair_sal,'b.'); axis([0 6 33.6 34.8]);
    title('LayerPair_sal');
    
end % if   

% -------------------


clear iprop Layer_mean IntLayer_prop binPair_prop IntLayerPair_prop
clear good_p good_surf PropName

disp('Calculating AreaPair ...')
deldistkm      = dobox_distance(lat,lon,'km');
distkm         = cumsum([0 deldistkm]);
Area_layer     = thick*diag(1000*deldistkm);
% remove next section elm 29.9.99
% note for clarification - elm 14/4/00 - area_tri is no longer addaed
% explicitly because thick is calculated from integrated frac   
%    if strcmp(reply,'y') == 1
%        add_tri = input('Add tri file:  ','s');
%        eval(['load ' add_tri]);
%         for iadd = 1:nstations -1
%           tt(iadd) = max(find(Area_layer(:,iadd) > 0));
%           Area_layer(tt(iadd),iadd) = ...
%                Area_layer(tt(iadd),iadd) + area_tri(iadd);
%         end %(for)
%     end %(if)
disp('Calculating Area_properties for')
disp('     mass')

Area_mass = Area_layer;


for iprop = 2:nproperties
   PropName = [charword( properties(iprop,:) )];
   disp(['     ' PropName])
   LayerPair_prop = ['LayerPair_'   PropName];

%  Area_prop  = Area_layer .* LayerPair_prop;
      VariableName = ['Area_' PropName];
      Expression   = ['Area_layer .* ' LayerPair_prop];
      command = [VariableName ' = ' Expression ';'];
      eval(command)
end %for   
clear iprop VariableName Expression command LayerPair_prop PropName
clear tt add_tri
clear TRUE FALSE

prompt_user = 0;  % do not prompt user to overwrite if sectfile passed
if exist('sectfile')==0
   sectfile = input('What filename to save this section? ','s');
   prompt_user = 1;
end %if

[sectfile,file_exists] = validfile(sectfile,'.mat');
if file_exists
   if prompt_user
       reply = input('File exists, overwrite (y/n) ? [n] ','s');
       if ~strcmp(lower(reply),'y')
           error('ABORT. No overwriting of existing file')
       end %if
   end %if
else
   % ok - batch mode
end %if
clear file_exists prompt_user


disp(['Saving file {' sectfile '} '])
command = ['save  ' sectfile];
setstr(command);
eval(command);

return

