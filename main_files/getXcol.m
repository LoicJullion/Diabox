%==========================================================================
% The diapycnal velocities have been added in with the get_diapvel function
% Now we refine the additional terms of the inverse model.
% Note that this function should be updated if one want to new terms to the
% inversion.
% As it stands, this routine:
%   - Adjusts the diapycnal velocities for some bugs
%   - Take into account the advection across isopycnals for properties
%   other than mass (PT, S).
%   - Load the air-sea fluxes and Ekman transport calculated by the 
%   makeairsea.m file 
%
% AUTHORS: Loic Jullion based on Rick Lumpkin's code.
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
%% Adjust diapycnal velocities
% Make sure w* terms for deepest 2 layers are consistent
% (bug identified 26 March 2004 - present in "readyforgauss"
% verion of A)
for propnum=1:nproperties;
    for boxnum=1:nboxes;
        A(indb(nlayers,propnum,boxnum),indw(nlayers-1,propnum,boxnum))= ...
        -A(indb(nlayers-1,propnum,boxnum), indw(nlayers-1,propnum,boxnum));
    end;
end;
% This could be fixed in Xcols.mat when I recalculate these terms
%--------------------------------------------------------------------------

disp('Conserved properties are hard-coded. ');
disp('  CHECK: property 1 should be mass,');
disp(['                       and is ',properties(1,:)]);
disp('  CHECK: property 2 should be salam,');
disp(['                       and is ',properties(2,:)]);
disp('  CHECK: property 3 should be potmp or temp,');
disp(['                       and is ',properties(3,:)]);
meanprop_nonans=meanprop;
level_area_nonans=level_area;
for boxnum=1:nboxes;
  for propnum=1:nproperties;
      meanprop_nonans(:,propnum,boxnum) = fillnans(meanprop(:,propnum,boxnum),0);
  end;
end;

q                    = find(isnan(level_area));
level_area_nonans(q) = 0;

% Include adjustments to density conservation (i.e. w*_rho terms) in 
% budgets of other properties - this represents advection across a mean 
% property gradient driven by net buoyancy fluxes, in addition to the 
% diffusive (eddy) fluxes parameterized by the w* terms.
for boxnum=1:nboxes;
  for propnum=2:nproperties;
    for laynum=1:nlayers;
      A(indb(laynum,propnum,boxnum),indw(:,1,boxnum))= ...
	  A(indb(laynum,propnum,boxnum),indw(:,propnum,boxnum));
    end
  end
end

%% Include Air sea fluxes 
%  1) Define matrices which contain the surface fluxes and values
%
%  AirSea(layer,property,box): air-sea input
%  eAirSea(layer,property,box): error bar range of AirSea
%  EkmanI(layer,property,box): internal (within each box) Ekman flux
%  EkmanE(layer,property,box,section): external (across-section) Ekman flux
%  SSValue(layer,property,box): sea surface area, temp, salinity
%  Si_Input(box): river input of silica (kmol/s)

Fdens       = NaN*ones(nlayers+1,nboxes);
AirSea      = NaN*ones(nlayers,nproperties,nboxes);
eAirSea     = AirSea;EkmanI=AirSea; SSValue=AirSea;
EkmanE      = NaN*ones(nlayers,nproperties,nboxes,nsections);
EkmanAdjust = zeros(nsections,1);
Si_Input    = NaN*ones(nboxes,1);

load dir_loc.mat as_dir sect_dir
for boxnum=1:nboxes;
  eval(['load ' as_dir 'airsea',num2str(boxnum),' M* F* Ek*' ...
	' gamma2 SS* si_flux']);

  Fdens(:,boxnum) = Fm(1,:)';

  AirSea(:,1,boxnum)   = Ffw(1,:);
  AirSea(:,2,boxnum)   = zeros(size(Ffw(1,:)));
  AirSea(:,3,boxnum)   = Fh(1,:);

  eAirSea(:,1,boxnum)  = Ffw(2,:);
  eAirSea(:,2,boxnum)  = Mm(2,:);
  eAirSea(:,3,boxnum)  = Fh(2,:);

  EkmanI(:,1,boxnum)   = Mmek;
  EkmanI(:,2,boxnum)   = Msek;
  EkmanI(:,3,boxnum)   = Mhek;

  EkmanE(:,1,boxnum,:) = reshape(EkmanInM, [nlayers,1,1,nsections]);
  EkmanE(:,2,boxnum,:) = reshape(EkmanInS, [nlayers,1,1,nsections]);
  EkmanE(:,3,boxnum,:) = reshape(EkmanInH, [nlayers,1,1,nsections]);

%   EkmanAdjust = max([EkmanAdjust, ...
%                      EkmanInM_range./abs(nansum(EkmanInM)'+1) ...
%                      ]')';

  SSValue(:,1,boxnum) = SSArea;
  SSValue(:,2,boxnum) = SSSbar;
  SSValue(:,3,boxnum) = SSTbar;

  Si_Input(boxnum)    = si_flux;

end;
clear Mm* Ffw* Fh* Ms* Mh* EkmanInM EkmanInS EkmanInH;
clear SSTbar SSSbar SSArea;

% Find interior sections (shared by two boxes)
% and insure net Ekman flux out of one box equals flux
% into the other.   This won't be true by default, due to
% extrapolation of winds and changes in surface gamma
% across the section.  
% THIS ROUTINE ASSUMES A SECTION CAN BE SHARED BY TWO AND ONLY TWO BOXES.
sectint=find(sum(geometry,1)==0);  % identify shared sections
for i = 1:length(sectint);
  EkM  =[];EkS=[];EkH=[]; % matrices to store means
  for boxnum = 1:nboxes;  % calculate mean Ekman flux
      EkM = [EkM,EkmanE(:,1,boxnum,sectint(i)).* ...
             geometry(boxnum,sectint(i))];
      EkS = [EkS,EkmanE(:,2,boxnum,sectint(i)).* ...
             geometry(boxnum,sectint(i))];
      EkH = [EkH,EkmanE(:,3,boxnum,sectint(i)).* ...
             geometry(boxnum,sectint(i))];
  end;
  % Resulting EkM(nlayers x nboxes) is the Ekman 
  % transport in each layer, into each box, for section
  % sectint(i).  All but two of the boxes are zero;
  % the two which share sect(sectint) now have transports
  % with the same sign, because they have been multiplied
  % by the corresponding value in "geometry" (one is +1,
  % the other is -1).  Sum the values over the box index
  % and divide by two to get the mean Ekman transport
  % into/out of the bounding boxes.
  EkM = sum(EkM,2)/2;
  EkS = sum(EkS,2)/2;
  EkH = sum(EkH,2)/2;
  for boxnum=1:nboxes;		% assign mean to each section
    if (geometry(boxnum,sectint(i)))  % section borders box
      EkmanE(:,1,boxnum,sectint(i))= ...
	EkM/geometry(boxnum,sectint(i));
      EkmanE(:,2,boxnum,sectint(i))= ...
        EkS/geometry(boxnum,sectint(i));
      EkmanE(:,3,boxnum,sectint(i))= ...
        EkH/geometry(boxnum,sectint(i));
    end;
  end;
end;
clear sectint EkM EkS EkH;

% Identify which layers exist in each box using the bounding hydrography 
% around each box and the interior air-sea forcing.
defined_layers = zeros(length(glayers),nboxes);
for boxnumber=1:nboxes;
  q=find(abs(geometry(boxnumber,:)));
  for i=1:length(q);
      eval(['load ',sectfiles(q(i),:),' thick']);
      thick=sum(thick')';
      q2=find(thick>0 | (AirSea(:,3,boxnumber)~=0));
      defined_layers(q2,boxnumber)=1;
  end; % loop over sections bounding box
end; % loop over boxnumber
clear thick q q2;


%% 2) Define the a priori diapycnal advection 
% 
% if ~includemixing_ap
   drD_ap=zeros(nlayers+1,nproperties,nboxes);
   edrD_ap=drD_ap;
% end;
drD_apDivA=drD_ap;
edrD_apDivA=edrD_ap;
% convert drD_apDivA to units of m/s (same as w* terms)
for boxnum=1:nboxes;
  for propnum=1:nproperties;
    divfact=meanprop_nonans(:,propnum,boxnum).* ...
	level_area(:,boxnum);
    drD_apDivA(2:nlayers,propnum,boxnum)= ...
        drD_apDivA(2:nlayers,propnum,boxnum)./divfact;
    edrD_apDivA(2:nlayers,propnum,boxnum)= ...
     	edrD_apDivA(2:nlayers,propnum,boxnum)./abs(divfact);
  end;
end;
q=find(isnan(drD_apDivA));
drD_apDivA(q)=0;
edrD_apDivA(q)=0;
% define W_ap: given by prior Fdens and mixing
W_ap=zeros(nlayers+1,nboxes);
for boxnum=1:nboxes;
  ind=2:nlayers;
  W_ap(ind,boxnum)=-Fdens(ind,boxnum)./ ...
	level_area(:,boxnum) + ...
	drD_apDivA(ind,1,boxnum);
end;
q=find(isnan(W_ap));
W_ap(q)=0;
clear drD_apDivA;
