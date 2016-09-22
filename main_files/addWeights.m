%==========================================================================
% The a priori uncertainties on:
%   - the conservation equations (row)
%   - the unknowns variables (columns)
% of the model are added in this mfile. This file needs to be updated by
% the user in order to reflect the importance of the different terms of the
% inverse model.
% The file needs also to be updated if new terms are added to the inversion
%
% OUTPUT:
%   - Rxx: Covariance matrix of the unkowns terms.
%   - Rnn: Model Covariance matrix.
% Both are usually taken as diagonal which means that all the equations
% (unknowns) are statistically independant from each other. In other word,
% varibility of one term does not affect the variability of the other
% terms. We know this is a poor approximation but based on in situ data, it
% is very complicated to obtain an estimation of the non-diagonal terms.
% See Wunsch, 1996).
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
%% Setup a priori errors
masslayererror=1;	% layer mass error (Sv)
ms_mf=2;      		% minsalamerror (psu); 35 in orig.
mindiaflux=10;		% floor to allowed range of 
                    % dia.vol.flux (Sv) driven by mixing

includemixing_ap = 0;	% set to 1 to include drD_ap

SILCATerror=500;	% silicate conserved at this rate (kmol/s)
% original: 10.  Robbins and Toole: 100.  Ganachaud et al: 700

masslayererror=masslayererror*ones(nlayers,1); % One error per layer

masslayererror=masslayererror.*1e6; % transform in Sv

mindiaflux=mindiaflux*1e6; % transform in Sv

%% % 1) Errors on the conservation equations
% allocate memory for model error
% nlayers + 1 because we add a row for the full depth transport
modelerror = NaN*ones(nlayers+1,nproperties,nboxes);
load glevels.mat glevels
% Loop over all the boxes
for boxnum=1:nboxes
  % 1a) mass errors 
  % set layer errors to masslayererror
  modelerror(1:nlayers,1,boxnum) = masslayererror;

  % assign masstotalerror 
  modelerror(nlayers+1,1,boxnum)=masstotalerror(boxnum);
  %!!!!!!!!!!!!!!!!
  % There is a limitation here: The high errors are assigned to the first 
  % layers. This works well with only a few boxes where density layers in
  % each boxes are more or less the same. When doing a global inversion,
  % the first layers (very light), will only be present in the tropical
  % boxes. This means that for mid-latitude and subpolar boxes, even the
  % first few layers will have relatively low uncertainties.
  %!!!!!!!!!!!!!!!!

  % identify the sections surrounding this box
  sects = find(geometry(boxnum,:));
  % loop through sections, loading layer-mean properties
  potmp_box_layer=[];salam_box_layer=[];
  for sectnum = 1:length(sects);
        eval(['load ',sectfiles(sects(sectnum),:), ...
	'   LayerPair_potmp LayerPair_salam;']);
      potmp_box_layer = [potmp_box_layer,LayerPair_potmp];
      salam_box_layer = [salam_box_layer,LayerPair_salam];
  end;
  % remove 0s (empty pairs in LayerPair_*)
  % and calculate layer std.deviations.
  q = find(~potmp_box_layer);potmp_box_layer(q)=NaN;
  q = find(~salam_box_layer);salam_box_layer(q)=NaN;
  
  mean_salam = abs(nanmean(salam_box_layer')');
  mean_potmp = abs(nanmean(potmp_box_layer')');
  std_potmp  = fillnans(nanstd(potmp_box_layer')',0);
  std_salam  = fillnans(nanstd(salam_box_layer')',0);

  % assign (abs(mean) + 2 std)*mass flux error to modelerror
  modelerror(1:nlayers,2,boxnum)= (abs(mean_salam)+2*std_salam).* ...
	modelerror(1:nlayers,1,boxnum);
  % set floor to error
  q = find(modelerror(1:nlayers,2,boxnum)< ms_mf/1e3*modelerror(1:nlayers,1,boxnum));
  modelerror(q,2,boxnum) = ms_mf/1e3*modelerror(q,1,boxnum);

  % vert. integrated error: max layer error times
  % ratio of total mass error to layer mass error
  %q=find(modelerror(1:nlayers,2,boxnum)== ...
        %max(modelerror(1:nlayers,2,boxnum)));
  %modelerror(nlayers+1,2,boxnum)= ...
        %modelerror(q(1),2,boxnum) ...
        %* modelerror(nlayers+1,1,boxnum) ...
        %/ modelerror(q(1),1,boxnum);
  % new method: sum the layer errors as if they
  % were independent, and scale by
  % ratio of total mass error to layer mass errors
  %modelerror(nlayers+1,2,boxnum)= ...
  %      sqrt(sum(modelerror(1:nlayers,2,boxnum).^2))  ...
  %      *masstotalerror(boxnum)/sqrt(sum( ...
  %     modelerror(1:nlayers,1,boxnum).^2));
  % newest method: total salt error is dominated by
  % upper layers, so overly-constrained net is
  % equivalent to overly-constrained upper layers.
  % This is a weaker constraint than the "new method" 
  % described above.
  modelerror(nlayers+1,2,boxnum)= ...
      sqrt(sum( (modelerror(1:nlayers,2,boxnum) - ...
      min(modelerror(1:nlayers,2,boxnum))) .^2 )) + min(modelerror(1:nlayers,2,boxnum));

  % if total salam error is < ms_mf/1e3*masstotalerror,
  % set it to this value
  if (modelerror(nlayers+1,2,boxnum) < ms_mf/1e3*masstotalerror(boxnum))
    modelerror(nlayers+1,2,boxnum) = ms_mf/1e3*masstotalerror(boxnum);
  end;

  %%%% potential temperature errors %%%%
  % error in PW: meanprop(:,3,boxnum)/(1e15/3990/1026).
  % assign (abs(mean) + 2 std)*mass flux error to modelerror
  modelerror(1:nlayers,3,boxnum)= ...
        (abs(mean_potmp)+2*std_potmp) .* modelerror(1:nlayers,1,boxnum);

  % vert. integrated error: max layer error times
  % ratio of total mass error to layer mass error
  q=find(modelerror(1:nlayers,3,boxnum)== max(modelerror(1:nlayers,3,boxnum)));
  modelerror(nlayers+1,3,boxnum)= modelerror(q(1),3,boxnum) ...
	* modelerror(nlayers+1,1,boxnum) / modelerror(q(1),1,boxnum);

  % new method: sum the layer errors as if they
  % were independent, and scale by
  % ratio of total mass error to layer mass errors
  %modelerror(nlayers+1,3,boxnum)= ...
  %	sqrt(sum(modelerror(1:nlayers,3,boxnum).^2))  ...
  %	*masstotalerror(boxnum)/sqrt(sum( ...
  %	modelerror(1:nlayers,1,boxnum).^2));
  % newest method: total heat error is dominated by
  % upper layers, so overly-constrained net is
  % equivalent to overly-constrained upper layers.
  % This is a weaker constraint than the "new method" 
  % described above.
  modelerror(nlayers+1,3,boxnum)= ...
      sqrt(sum( (modelerror(1:nlayers,3,boxnum) - ...
      min(modelerror(1:nlayers,3,boxnum))) .^2 )) + min(modelerror(1:nlayers,3,boxnum));
  
  
  % January 2007: very tight constraint on net heat; this
  % can be accomodated by adjusting Ekman and/or heat fluxes.
  %modelerror(nlayers+1,3,boxnum)=min( ...
  %    [modelerror(nlayers+1,3,boxnum) .05*(1e15/3990/1026)]);

  eval(['clear *_atl',num2str(boxnum)]);
  clear mean_salam* std_salam mean_potmp* std_potmp;

end;

% allocate memory for the noise and model
% variance matrices
Rnn=eye(mA,mA);
Rxx=eye(nA,nA);

% assign squared values to error matrix Rnn
for i=1:nlayers+1;
  for j=1:nproperties;
    for k=1:nboxes;
      Rnn(indb(i,j,k),indb(i,j,k))=modelerror(i,j,k)^2;
    end;
  end;
end;


%% 2) Set a priori values for unknowns %

%   2a) a priori reference velocities
Rxx(indv,indv)=diag(vmag.^2);
% 
% if (0)
%   % find where station separation distance
%   % is less than corrdist (km); scale the cross-
%   % correlation elements of Rxx for these
%   % pairs by the a priori vmag values,
%   % scaled by ( 1-dist/corrdist).
%   corrdist=20;
%   for stnpair=1:length(lon);
%     dist=NaN*ones(size(lon));
%     for j=1:length(lon);
%       dist(j)=sw_dist([lat(stnpair) lat(j)], ...
% 	[lon(stnpair) lon(j)],'km');
%     end;
%     q=find(dist>0 & dist<corrdist);
%     for j=1:length(q);
%       Rxx(indv(q(j)),indv(stnpair))= ...
% 	vmag(q(j))*vmag(stnpair) * (1-dist(q(j))/corrdist);
%       Rxx(indv(stnpair),indv(q(j)))= ...
% 	vmag(q(j))*vmag(stnpair) * (1-dist(q(j))/corrdist);
%     end;
%   end;
% end;

%   2b) a priori diapycnic velocities

for boxnum=1:nboxes;
    
  wmag = NaN*ones(nlayers-1,nproperties);
  for propnum = 1:nproperties

    pfluxrel = -sum(D(indb(1:nlayers,propnum,boxnum),indv),2);
    pfluxek  = sum(EkmanE(1:nlayers,propnum,boxnum,:),4);
    Psi      = [flipud(cumsum(flipud(pfluxrel+pfluxek)));0];
    % include direct air-sea input
    Psi = Psi+[flipud(cumsum(flipud(AirSea(:,propnum,boxnum))));0];
    % include a priori diapycnal advection
    % (note: "Psi" becomes Psi+drD-F at this step)
    wap=[W_ap(:,boxnum), W_ap(:,boxnum)];
    if (propnum>1);  % include adjustments to advection
      wap=wap+[wadj,-wadj];
    end;
    %note that -W_ap(:,1).*[0;level_area_nonans(:,1);0]=Fdens(:,1);
    Psi = [Psi,Psi]-[[0,0];wap(2:nlayers,:).* ...
	(level_area_nonans(:,boxnum)*ones(1,2) ).* ...
	(meanprop_nonans(:,propnum,boxnum)*ones(1,2));[0,0]];
    clear wap;
    
    %----------------------------------------------------------------------
    % From Rick's code
    % adjusted Psi will conserve net property, i.e. go to zero
    % at gamma2(1).  Assume these adjustments are evenly-
    % distributed across the layers
    % NOTE: if this is not done for mass, then accumulated 
    %  errors in the relative+ekman can produce large positive
    %  values of Psi in the lightest layers - this is reduced
    %  by the large negative (to lighter layers) F_as, producing
    %  small values of wmag when large ones are needed to
    %  counter F_as (Niiler-Stevenson balance)
    % Note by Loic: Not sure I really understand this bit of the code. It
    % may matter for a global inversion but probably not for regional
    % inversions
    if (0)
      mPsi=nanmean(Psi')';
      q=min(find(mPsi==0 & gamma2>27));
      q2=max(find(~abs(diff(mPsi)) & glayers<27))+1;
      if ~length(q2) q2=1; end;
      mPsi=fillnans((interp1([q2 q],[Psi(q2) Psi(q)], ...
        1:length(gamma2)))',0) *ones(1,2);
      Psi=Psi-mPsi;
      clear mPsi;
    end
    %----------------------------------------------------------------------
    
    if (propnum==1) 
      wadj    = Psi(:,1)./[NaN;level_area(:,boxnum);NaN];
      q       = find(isnan(wadj));
      wadj(q) = 0;
    end;
  
    diafluxmag = abs(-Psi);
    diafluxmag = nanmax(diafluxmag')';

    % drop endpoints to get glevels grid
    diafluxmag      = diafluxmag(2:length(diafluxmag)-1);
    wmag(:,propnum) = diafluxmag./level_area(:,boxnum) ./abs(meanprop(:,propnum,boxnum));
  end % loop over propnum

  % Allow for a little "slop" in the system: don't force
  % each and every level to have the minimum mixing according
  % to the a priori choices, so that the system has more
  % freedom to adjust from those a priori values and find
  % a better global minimum for mixing.  This "slop" is
  % mindiaflux (minimum diapycnal volume flux) for mass,
  % mindiaflux*(magprop) for other properties.
  
  %% experimental code 14 Nov 2003: increase this slop
  % in outcropping layers (identified by nonzero heat flux).
  % Allowing large mixing in these layers is consistent
  % with strong lateral eddy fluxes across isopycnals
  % in the mixed layer.
  if (1)
      q = max(find(AirSea(:,3,boxnum)~=0));
      if (q>nlayers-1) q=nlayers-1; end;
      mindiafluxL = mindiaflux*ones(nlayers-1,1);
      mindiafluxL(1:q) = 2*mindiafluxL(1:q);
      sarea            = level_area(:,boxnum) * ones(1,nproperties) ...
                             .*abs(meanprop(:,:,boxnum));
      magprop = nanmean(wmag.*sarea);
      magprop = magprop./magprop(1);
      wmag    = wmag+(mindiafluxL*magprop)./sarea;
      clear mindiafluxL q;
  end
    
    wmag(find(isnan(wmag))) = 1e-99;
    

  % include error bar on a priori mixing
  wmag = wmag+edrD_apDivA(2:nlayers,:,boxnum);

  % assign wmag^2 to Rxx
  for propnum = 1:nproperties
    Rxx(indw(:,propnum,boxnum),indw(:,propnum,boxnum)) = diag( wmag(:,propnum).^2 );
  end

  % air-sea adjustment terms
  % heat
  tvals = NaN*ones(size(glevels));
  for ti = 1:length(glevels)
    tvals(ti) = max([eAirSea(ti,3,boxnum) eAirSea(ti+1,3,boxnum)].^2);
  end
  q = find(tvals>0);q2=find(tvals==0);
  tvals(q2) = .01*min(tvals(q));
  Rxx(indF(:,1,boxnum),indF(:,1,boxnum)) = diag( tvals );

  % freshwater
  tvals = NaN*ones(size(glevels));
  for ti = 1:length(glevels)
    tvals(ti) = max([eAirSea(ti,1,boxnum) eAirSea(ti+1,1,boxnum)].^2);
  end
  q = find(tvals>0);q2=find(tvals==0);
  tvals(q2) = .01*min(tvals(q));
  Rxx(indF(:,2,boxnum),indF(:,2,boxnum)) = diag( tvals );

end	% looping over boxnum

% a priori adjustments to net Ekman transport
if (length(EkmanAdjust)==1)
  Rxx(indE,indE)=diag( EkmanAdjust * ones(size(indE))).^2;
else
  Rxx(indE,indE)=diag( EkmanAdjust ).^2;
end;
