%==========================================================================
% gauss.m
%
% This is the main routine which is going to put all the terms of the
% inversion calculated in the previous steps together and run the inverion
% Run after correcting y for air-sea forcing and using DOBOX to add 
% diapycnal flux terms to A, x.
% 
% The original gauss.m routine contained all the elements in one big file.
% In the new version, the different sections have been placed in separate
% mfiles. This should make it easier to undertand.
%
% Use Gauss-Markov estimation to solve the
% system Ax + n = y.
%
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
%% Define some variables
% fluxunits: units of terms in b, A*x (after dividing
%	by fluxdivfact).  Untouched units are
%       SI: volume is  salam is ,
%	potmp is dgC m^3/s.  Use a mean value of 
%	 for potemp to get PW.
fluxunits   = ['   Sv  '; ...
               ' psu Sv'; ...
               '   PW  ']; 
% fluxdivfact allows to convert the inversion results into the units fluxunits           
fluxdivfact = [1e6; ... % m^3/s,
               1e-3*1e6; ... % 1e-3 psu m^3/s
               1e15/3985/1026]; % 1e15/(rho*C_p) (PW)
propnames   = ['mass';'salt';'temp'];

propunits   = ['10^{6} km^2'; ...
	           '    psu    '; ...
               '  ^\circC  '];
propdivfact = [1e12;1e-3;1];

propcolr    = [[0 0 1];[0 .5 0];[1 0 0]];

load dir_loc.mat sect_dir refvel_dir;
% load positions of station pairs
lon2=[];lat2=[];theta2=[];sectlabel=[];R2=[];vmag=[];
for i=1:nsections;
  eval(['load ' sectfiles(i,:),' lon lat']);
  sectt = [refvel_dir 'refvel_', sectfiles(i,length(sect_dir)+1: ...
            size(sectfiles,2))];
  eval(['load ',sectt,' Vmag']);

  q  = find(abs(diff(lon))>90);	% +-180
  l2 = (lon(1:length(lon)-1)+diff(lon)/2)';
  lon(find(lon<0)) = lon(find(lon<0))+360;
  l2(q) = (lon(q-1)+lon(q))/2;
  if (l2(q)>180) l2(q)=l2(q)-360; end; % Convert to -180 180
  lon2      = [lon2;l2]; clear l2;
  lat2      = [lat2;(lat(1:length(lat)-1)+diff(lat)/2)'];
  [R,theta] = sw_dist(lat,lon,'km');
  R2        = [R2;R'];
  theta2    = [theta2;theta'];
  sectlabel = [sectlabel;i*ones(length(lon)-1,1)];
  vmag      = [vmag;Vmag'];
end;	

lon   = lon2;clear lon2;
lat   = lat2;clear lat2;
R     = R2*1000;clear R2;	% R [m]
theta = theta2*2*pi/360; clear theta2; % theta [rad]	
			
clear Vmag sectt;

%% Define the indices for the variables.
% This is important since it allows to track the location of the different 
% terms. Particularly important for inversion with a large number of boxes. 
% For 1 or 2 boxes, it is possible to track the variables manually, but
% this quickly becomes problematic when the number of boxes increases.

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
% NOTE THAT THIS MFILE NEEDS TO BE UPDATED IN ONE WISHES TO ADD NEW TERMS.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

[indv,indw,indF,indE,indb,meanprop,level_area] = defindices(A,nAbase,nA,...
    nboxes,nproperties,nlayers,nsections);
%% Load extra variables (air-sea fluxes, Ekman...)
% In the first part of DIABOX, the diapycnal fluxes have already been
% included.
addXterm = 0;
if addXterm == 1
    % 1) First load the new terms that have been previously calculated in 
    % makeairsea.m 
    getXcol;

    % 2) Now we add the new terms loaded in 1) to the terms of the inversion. 
    updateA_b;

    % !!!!!! NOTE THAT THESE 2 MFILES NEED TO BE UPDATED IN ONE WISHES TO ADD 
    % NEW TERMS. !!!!!
else
    disp('No new extra unkowns have been added. Inversion is only using diapycnal fluxes');
end
%% Add uncertainties
% !!!!!! NOTE THAT THIS MFILE NEEDS TO BE UPDATED ACCORDING TO THE WEIGHTS
% THE USER WISHES TO ADD TO THE DIFFERENT TERMS OF THE INVERSION.
% This is an important par since the uncertainty will reinforce, decrease
% the relative importance of the different terms of the inverse model.
addWeights;
%% Add extrat constraints 
xconst = 0;
if xconst == 1
    add_xconstraints;
else
    disp('No extra constraints have been added')
end
%% Solve the system via Gauss-Markov estimation 
disp('Solving the set of equations.');
fprintf('  Inverting ... ');
time=cputime;
warning off MATLAB:nearlySingularMatrix
C=Rxx*A'/(A*Rxx*A'+Rnn);
warning on MATLAB:nearlySingularMatrix
fprintf('done in %5.0f cpu seconds.\n',cputime-time);
% alternative; takes longer in tests on lasalle
%E=A*Rxx*A'+Rnn; [L,U]=lu(E); C=Rxx*A'/U/L;

x=C*b;

% calculate noise, model uncertainty
n=b-A*x;
Pxx=Rxx-C*A*Rxx;
%Pnn=eye(size(Rnn))-A*C;Pnn=Pnn*Rnn*Pnn;
erx=sqrt(diag(Pxx));
clear C;

% define Ax.  A*x=sum(Ax,2).
Ax=A.*(ones(length(b),1)*x');
Aerx=A.*(ones(length(b),1)*erx');

Fdenscorr=NaN*ones(nlayers-1,nboxes);
eFdenscorr=Fdenscorr;
for k=1:nboxes;
  Fdenscorr(:,k)= Fcorr_multfact(:,1,k).*x(indF(:,1,k))+ ...
	Fcorr_multfact(:,2,k).*x(indF(:,2,k));
  eFdenscorr(:,k)=sqrt( ...
	(Fcorr_multfact(:,1,k).*erx(indF(:,1,k)) ).^2 + ...
        (Fcorr_multfact(:,2,k).*erx(indF(:,2,k)) ).^2 );
end;
Fdenscorr=[zeros(1,nboxes);Fdenscorr;zeros(1,nboxes)];
eFdenscorr=[zeros(1,nboxes);eFdenscorr;zeros(1,nboxes)];

% include E* terms in EkmanE
for k=1:nsections;
  EkmanE(:,:,:,k)=EkmanE(:,:,:,k)*(1+x(indE(k)));
end;

% calculate the final diapycnal advection (m/s)
W=W_ap;
for boxnum=1:nboxes;
  W(2:nlayers,boxnum)=W(2:nlayers,boxnum) ...
	-Fdenscorr(2:nlayers,boxnum)./ ...
		level_area(:,boxnum) ...
	+ x(indw(:,1,boxnum));
end;
q=find(isnan(W));W(q)=W_ap(q);

clear Fer Psi Psi2 RV XCW XRW Xbs Xcols Xrows 
clear Xzeros alpha arcPE beta binPair_vel
clear bl box_menu_opt box_menu_stop boxnum
clear col_beg col_end cons ct defum dg
clear diafluxmag directory eFer eqn_beg eqn_end
clear file_Xcols file_Xrows ind j jj k laynum
clear multfact nboxeqns nsurfs pfluxek
clear q q2 rho rho0 row_beg row_end rv
clear rvfilename sarea sbar sectname stdsalt
clear tbar ti tvals wmag wadj

load dir_loc.mat output_dir;
eval(['save ' output_dir 'gaussout.mat;']);

if ( max(abs(x)./sqrt(diag(Rxx)) )>1 | ...
	max(abs(n)./sqrt(diag(Rnn)) )>1)
  disp('*** inversion failed to satisfy a priori choices ***');
end;

%% Plot results 

figure(1);	% compare x and n to priors
compsoln2priors;

return;

figure(2)
propnum=1;
massnum=0;
plotlatflux;

figure(3);	% model solution velocities
clf;
propnum=1;
for boxnum=1:nboxes;
  subplot(2,nboxes,boxnum);
  plotdiaflux;
end;
subplot(212);
plotrefvel;

