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
% AUTHOR:  Phil Morgan 12-03-92 Modified by Rick Lumpkin and Loic Jullion
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
   disp(['Loading data file {' sectfile ',_raw.mat}'])
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
Press2Row         = [bin_press [1:length(bin_press)]' ];
[values,com_rows] = lastgood( eval(['binPair_' PropName2]),NaN);

pelp1=fliplr(Press2Row);
com_press = interp1(pelp1(:,1),pelp1(:,2),com_rows);
clear values PropName2 pelp1

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

Layer_press=NaN*ones(nsurfs+2,npairs);
Area_layer=NaN*ones(nsurfs+1,npairs);

disp('Calculating thickness, area of Layers ... ')
thick=NaN*ones(size(surf_press,1)+1,length(lon)-1);
fprintf('  counting to %d: ',size(thick,2));

for i=1:length(lon)-1;
  fprintf('%d ',i);

  % extract values of surf_press from both casts
  sp1=surf_press(:,i);
  sp2=surf_press(:,i+1);

  q=find(isfinite(sp1)); if (length(q))
  if sum(sp1(q)~=sort(sp1(q)))
    fprintf('  gamma not monotonic in cast %d.\n',i);
    sp1(q)=sort(sp1(q));
  end;end;
  q=find(isfinite(sp2)); if (length(q))
  if sum(sp2(q)~=sort(sp2(q)))
    fprintf('  gamma not monotonic in cast %d.\n',i+1);
    sp2(q)=sort(sp2(q));
  end;end;

  % Add top and bottom level for layers 1, nsurfs+1.
  % Add NaNs to top if upper layer doesn't exist
  % or is at 0, otherwise add 0.  Same with bottom.
  if (sp1(1)>0) sp1=[0;sp1];
  else sp1=[NaN;sp1]; end;
  if (sp2(1)>0) sp2=[0;sp2];
  else sp2=[NaN;sp2]; end;
  if (last(sp1)<bottom(i)) sp1=[sp1;bottom(i)];
  else sp1=[sp1;NaN]; end;
  if (last(sp2)<bottom(i+1)) sp2=[sp2;bottom(i+1)];
  else sp2=[sp2;NaN]; end;

  % In cases where only the uppermost or lowermost
  % layer exists in a cast, spN is all NaNs.  Fix.
  if (~sum(isfinite(sp1)) & sum(isfinite(sp2)) )
    if (mean(find(isfinite(sp2)))>((nsurfs+2)/2)) 
      sp1(nsurfs+1)=0;sp1(nsurfs+2)=bottom(i);
      fprintf('  assigning deepest layer to cast %d.\n',i);
    else
      sp1(1)=0;sp1(2)=bottom(i);
      fprintf('  assigning layer 1 to cast %d.\n',i);
    end;
  elseif (~sum(isfinite(sp2)) & sum(isfinite(sp1)) )
    if (mean(find(isfinite(sp1)))>((nsurfs+2)/2))
      sp2(nsurfs+1)=0;sp2(nsurfs+2)=bottom(i+1);
      fprintf('  assigning deepest layer to cast %d.\n',i+1);
    else
      sp2(1)=0;sp2(2)=bottom(i+1);
      fprintf('  assigning layer 1 to cast %d.\n',i+1);
    end;
  elseif (~sum(isfinite(sp1)) & ~sum(isfinite(sp2)) )
    disp('Cannot identify layers in dosection.m');
    keyboard;
  end;
  
  % should have min(diff(spN))>0.  Fix if not happening
  q=find(diff(sp1)==0);
  if (length(q))
      disp('WARNING!  diff(sp1) has a zero\n')
      if q<(nsurfs/2)
        sp1(q+1)=.5*(sp1(q)+sp1(q+2));
      else
        sp1(q)=.5*(sp1(q-1)+sp1(q+1));
      end;
  end;
  q=find(diff(sp2)==0);
  if (length(q))
      disp('WARNING! diff(sp2) has a zero - \n')
      if q<(nsurfs/2)
          sp2(q+1)=.5*(sp2(q)+sp2(q+2));
      else
        sp2(q)=.5*(sp2(q-1)+sp2(q+1));
      end;
  end;

  dx=sw_dist([lat(i) lat(i+1)], ...
        [lon(i) lon(i+1)],'km')*1e3;
  % calculate area of each layer between casts
  % using routine interpcasts.m
  [A,LayP]=interpcasts( ...
	sp1,sp2,[bottom(i) bottom(i+1)],dx,nsurfs);
  Area_layer(:,i)=A;
  thick(:,i)=A/dx;
  q=find(isnan(thick(:,i)));
  thick(q,i)=0;
  Layer_press(:,i)=LayP;
end;
fprintf('\n');
thick = change(thick,'==',NaN,0);
clear i A dx q lp lp2 sp1 sp2 LayP ans;

Layer_bin = NaN*ones(size(Layer_press));
for ipair = 1:npairs
   goodi=find( ~isnan(Layer_press(:,ipair)) );
   % length(goodi) = 0 (ok, no surfs))
   if length(goodi) > 1
      Layer_bin(goodi,ipair) = fix(interp1( ...
	Press2Row(:,1),Press2Row(:,2),Layer_press(goodi,ipair)) );
   elseif length(goodi) == 1
      error('dosection.m: must have more than 1 surface to define layers')
   end %if
end %for
clear Press2Row goodi

clear val1 val2 itop ibot

disp('Calculating mean property value in middle of LayerPair for ')

LayerPair_mass = ones(nsurfs+1,npairs);

for iprop = 2:nproperties
    PropName = [charword( properties(iprop,:) )];
    disp(['     ' PropName])
    binPair_property = eval(['binPair_'   PropName]);

    IntLayer_prop  = [];

    for ipair = 1:npairs
        IntLayerPair_prop   = NaN*ones(1,nlayers);
  
        % FIND VALID NEUTRAL SURFACES FOR A STATION PAIR. incl csurf and bottom
        good_p      = Layer_press(:,ipair); 
        good_surf   = find(~isnan(good_p)); 
  
        % INTEGRATE PROP DOWN A STATION PAIR BETWEEN VALID SURFACES
        if (~isempty(good_surf) & 0)
            F_good = [];
            F_good = integrate(bin_press,binPair_property(:,ipair),good_p(good_surf));
            % INSERT INTEGRATE FLUXES INTO VALID LAYERS, OTHERS LEFT AS NaN
            for igood = 1:length(F_good)
                IntLayerPair_prop(good_surf(igood))   = F_good(igood);
            end %for
        else 
            for igood=1:length(good_surf)-1;
                q=find(bin_press>=good_p(good_surf(igood)) & ...
                bin_press<=good_p(good_surf(igood+1)));
                if ~length(q) 
                    q2=find(bin_press>good_p(good_surf(igood+1)));
                    q=q2(1);
                end
                IntLayerPair_prop(good_surf(igood))=nanmean( ...
                binPair_property(q,ipair) );
            end
        end %if

        % APPEND INTEGRATION OF ONE STATION PAIR FLUXES INTO A SECTION OF FLUXES
        IntLayer_prop   = [IntLayer_prop IntLayerPair_prop'];

    end % loop over ipair
    clear F_good binPair_property igood good_p good_surf 
   
    Layer_mean = IntLayer_prop;
    Layer_mean = change(Layer_mean,'==',NaN,0);


    VariableName = ['LayerPair_' PropName];
    command = [VariableName ' = Layer_mean;'];
    eval(command)
   
end %for  

clear iprop Layer_mean IntLayer_prop binPair_prop IntLayerPair_prop
clear good_p good_surf PropName

deldistkm      = dobox_distance(lat,lon,'km');
distkm         = cumsum([0 deldistkm]);

disp('Calculating Area_properties for')
disp('     mass')
Area_mass = Area_layer;

for iprop = 2:nproperties
    PropName = [charword( properties(iprop,:) )];
    disp(['     ' PropName])
    LayerPair_prop = ['LayerPair_'   PropName];

    VariableName = ['Area_' PropName];
    Expression   = ['Area_layer .* ' LayerPair_prop];
    command = [VariableName ' = ' Expression ';'];
    eval(command)
end %for   
clear iprop VariableName Expression command LayerPair_prop PropName
clear tt add_tri
clear TRUE FALSE

clear q q2 

disp(['Saving file {' sectfile '} '])
command = ['save  ' sectfile];
setstr(command);
eval(command);

return

