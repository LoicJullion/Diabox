%==========================================================================
% Once all the extra terms have been succesfully loaded, the terms of the
% inversion can be updated:
%   - A: The main matrix containing all the terms of the equation
%   - b: The a priori estimate of the transport residual in each layer
%
% Note that this function should be updated if one want to new terms to the
% inversion.
% As it stands, this routine includes:
%   - the adjustment to the diapycnal velocities for some bugs
%   - the advection across isopycnals for properties other than mass (PT, S).
%   - the air-sea fluxes and Ekman transport calculated by the 
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
%% 3) Now include the air-sea formation, Ekman fluxes and a priori mixing in vector b
includeEkmanI=0;
disp('Including air-sea fluxes in vector b.');
for propnum=1:nproperties;
  for boxnum=1:nboxes;

    % include Ekman fluxes
    b(indb(1:nlayers,propnum,boxnum))= ...
        b(indb(1:nlayers,propnum,boxnum)) - ( ...
	includeEkmanI*EkmanI(:,propnum,boxnum)+ ...
        sum(EkmanE(:,propnum,boxnum,:),4)   );
    b(indb(nlayers+1,propnum,boxnum))= ...
        b(indb(nlayers+1,propnum,boxnum)) - nansum( ...
        includeEkmanI*EkmanI(:,propnum,boxnum)+ ...
        sum(EkmanE(:,propnum,boxnum,:),4)   );

    % include property convergence/divergence driven by a priori diapycnal 
    % advection
    b(indb(1:nlayers,propnum,boxnum))= ...
	b(indb(1:nlayers,propnum,boxnum)) ...
	-diff(        W_ap(:,boxnum).* ...
           [0;level_area_nonans(:,boxnum);0].* ...
           [0;meanprop_nonans(:,propnum,boxnum);0] );

    % include direct air-sea inputs (AirSea for salt is zero).
    b(indb(1:nlayers,propnum,boxnum))= ...
	b(indb(1:nlayers,propnum,boxnum)) - AirSea(:,propnum,boxnum);

    b(indb(nlayers+1,propnum,boxnum))= ...
	b(indb(nlayers+1,propnum,boxnum)) - nansum(AirSea(:,propnum,boxnum));

  end;
end;

%% 4) Include adjustments to air-sea forcing in the set of conservation equations
disp('Including adjustments to air-sea fluxes in A and x.');

A=[A,zeros(mA,(nlayers-1)*2*nboxes)];

Fcorr_multfact=NaN*ones(nlayers-1,2,nboxes);
% 

for k=1:nboxes;

  % add heat flux adjustments to heat budget
  % (remember, adjustments are defined on glevels
  %  grid; add 1/2 from level i-1, 1/2 from level i.)
  %  for top/bottom layers, use 1/2 value at closest level.)
  A(indb(1,3,k),indF(1,1,k))                     = .5;
  A(indb(2:nlayers-1,3,k),indF(1:nlayers-2,1,k)) = .5*diag(ones(nlayers-2,1));
  A(indb(2:nlayers-1,3,k),indF(2:nlayers-1,1,k)) = .5*diag(ones(nlayers-2,1));
  A(indb(nlayers,3,k),indF(nlayers-1,1,k))       = .5;
  
  % add sum of layer adjustments to total heat budget
  A(indb(nlayers+1,3,k),indF(:,1,k)) = ones(1,nlayers-1);

  % add freshwater flux adjustments to volume budget
  A(indb(1,1,k),indF(1,2,k))                     = .5;
  A(indb(2:nlayers-1,1,k),indF(1:nlayers-2,2,k)) = .5*diag(ones(nlayers-2,1));
  A(indb(2:nlayers-1,1,k),indF(2:nlayers-1,2,k)) = .5*diag(ones(nlayers-2,1));
  A(indb(nlayers,1,k),indF(nlayers-1,2,k))       = .5;
  
  % add sum of layer adjustments to total volume budget
  A(indb(nlayers+1,1,k),indF(:,2,k)) = ones(1,nlayers-1);

  %------------------------------------------------------------------------
  % Construct a density flux adjustment from the heat and freshwater 
  % flux adjustments, and include this advection in the property budgets.
  %
  % NOTE by LOIC: Not sure how significant this part of the code is. Not
  % used for now
  include_Fcorr_multfact = 0;
  if include_Fcorr_multfact == 1
      dg=diff(glayers);
      
      Tbar  = interp1(glayers,fillnans(SSValue(:,3,k),0), glevels,'linear');
      Sbar  = interp1(glayers,fillnans(SSValue(:,2,k),0), glevels,'linear');
      rho   = sw_dens(Sbar,Tbar,0);
      rho0  = sw_dens(zeros(size(Tbar)),Tbar,0);
      alpha = sw_alpha(Sbar,Tbar,0);
      beta  = sw_beta(Sbar,Tbar,0);

      % density flux from heat flux
      Fcorr_multfact(:,1,k)=-alpha.*rho./dg;

      % density flux from freshwater flux:
      % multiply by -1 (fw flux has opposite
      % sign of E-P flux) and don't forget
      % factor of SSS/(1-SSS/1000)
      Fcorr_multfact(:,2,k)=-beta.*rho./dg .*Sbar./(1-Sbar/1000);

  else
      Fcorr_multfact = zeros(nlayers-1,2,nboxes);

      % Loop through layers.  For each layer, add the
      % density flux into the layer at level i-1 and
      % remove the flux out of the layer at level i.
      for i=1:nlayers-1	% remove flux at glevel(i)
        for propnum=1:nproperties
            A(indb(i,propnum,k),indF(i,1,k))= A(indb(i,propnum,k),indF(i,1,k)) ...
              -Fcorr_multfact(i,1,k).* meanprop_nonans(i,propnum,k);
            A(indb(i,propnum,k),indF(i,2,k))= A(indb(i,propnum,k),indF(i,2,k)) ...
              -Fcorr_multfact(i,2,k).* meanprop_nonans(i,propnum,k);
        end
      end
      
      for i=2:nlayers	% add flux at glevel(i-1)
        for propnum=1:nproperties;
            A(indb(i,propnum,k),indF(i-1,1,k))= A(indb(i,propnum,k),indF(i-1,1,k)) ...
              + Fcorr_multfact(i-1,1,k).* meanprop_nonans(i-1,propnum,k);
          A(indb(i,propnum,k),indF(i-1,2,k))= A(indb(i,propnum,k),indF(i-1,2,k)) ...
            + Fcorr_multfact(i-1,2,k).* meanprop_nonans(i-1,propnum,k);
        end
      end
      % by construction, these terms sum to zero for
      % vertically-integrated property balances.

      % done adding density flux adjustments
  end % if
  %------------------------------------------------------------------------
end;
nA=size(A,2);

%% 5) Include adjustments to air-sea forcing in the set of conservation equations 

disp('Including adjustments to Ekman transports in A and x.');
A=[A,zeros(mA,nsections)];
for sectnum=1:nsections
    for i=1:nlayers+1
        for j=1:nproperties
            for k=1:nboxes
                if (i<nlayers+1)
                    A(indb(i,j,k),indE(sectnum))= EkmanE(i,j,k,sectnum);
                else
                    A(indb(i,j,k),indE(sectnum))= nansum(EkmanE(:,j,k,sectnum));
                end
            end
        end
    end
end

nA=size(A,2);
disp(' ');
