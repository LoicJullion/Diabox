function [drD_ap,edrD_ap]=makemixing_ap;

% function [drD_ap,edrD_ap]=makemixing_ap
%
% Transformation equation: A=F-drD.
% Formation equation: A=-Psi.
%
% At very light densities, A -> 0 and F is balanced
% by mixing (Niiler and Stevenson layers).
% An inverse model should include this necessary
% mixing in the a priori state.
%
% This routine determines the a priori mixing by
% plotting the formation equation terms for a box,
% highlighting where F+Psi=0.
% The user clicks on one of these
% points, and at all lighter densities an a priori
% mixing drD_ap is given by (F+Psi).  The a priori
% mixing at greater densities is zero, i.e. the
% inverse model must adjust Psi, F and drD to
% achieve property conservation.
%
% This a priori balance implies an a priori
% net diapcynal transformation WA (defined such that
% WA>0 implies transformation towards lighter layers).
% WA should act on the other properties according to 
% the divergence of WA*Cbar, Cbar=mean value of
% property C on the isopycnals.
%
% A priori budgets in each layer:
% 	(pfluxrel+pfluxek)-diff(F-drD_ap)
%
% If property j is balanced exactly by the a priori
% transports pfluxrel+pfluxek,
% then pfluxrel+pfluxek+diff(drD_ap(:,i,boxnum))=0.

clf;

clear all;load readyforgauss;
defindices;
EkmanE=NaN*ones(nlayers,nproperties,nboxes,nsections);
EkmanI=NaN*ones(nlayers,nproperties,nboxes);
AirSea=NaN*ones(nlayers,nproperties,nboxes);
drD_ap=zeros(nlayers+1,nproperties,nboxes);
WA=zeros(nlayers+1,nboxes);

load dir_loc.mat as_dir
for boxnum=1:nboxes;
for propnum=1:nproperties;

eval(['load ' as_dir 'airsea',num2str(boxnum),' M* F* Ek* gamma2 SS*']);
EkmanE(:,1,boxnum,:)=reshape(EkmanInM, ...
        [nlayers,1,1,nsections]);
EkmanE(:,2,boxnum,:)=reshape(EkmanInS, ...
        [nlayers,1,1,nsections]);
EkmanE(:,3,boxnum,:)=reshape(EkmanInH, ...
        [nlayers,1,1,nsections]);

EkmanI(:,1,boxnum)=Mmek;
EkmanI(:,2,boxnum)=Msek;
EkmanI(:,3,boxnum)=Mhek;

AirSea(:,1,boxnum)=Mm(1,:);
AirSea(:,2,boxnum)=Ffw(1,:);
AirSea(:,3,boxnum)=Fh(1,:);

Spbar=fillnans([meanprop(1,2,boxnum);meanprop(:,2,boxnum); ...
	meanprop(nlayers-1,2,boxnum)],0);
Tbar=fillnans([meanprop(1,3,boxnum);meanprop(:,3,boxnum); ...
	meanprop(nlayers-1,3,boxnum)],0);

pfluxrel=-sum(D(indb(1:nlayers,propnum,boxnum),indv),2);
pfluxek=sum(EkmanE(1:nlayers,propnum,boxnum,:),4);
Psi=[flipud(cumsum(flipud(pfluxrel+pfluxek)));0];
Aek=[flipud(cumsum(flipud( ...
	EkmanI(1:nlayers,propnum,boxnum))));0];

if (propnum==1)
  % adjust Psi to go to zero (conserve net volume)
  q=min(find(Psi==0 & gamma2>27));
  q2=max(find(~abs(diff(Psi)) & glayers<27))+1;
  if ~length(q2) q2=1; end;
  Psi=Psi-fillnans((interp1([q2 q],[Psi(q2) Psi(q)], ...
        1:length(gamma2)))',0);
end;
if (propnum==1)		% mass: buoy.flux convergence + fw flux
  F=[flipud(cumsum(flipud(AirSea(:,1,boxnum)))) + ...
	flipud(cumsum(flipud(AirSea(:,2,boxnum))));0];
elseif (propnum==2);	% salt: no air-sea input
  F=zeros(size(Psi));
else			% heat
  F=[flipud(cumsum(flipud(AirSea(:,3,boxnum))));0];
end;

drD=F+Psi;

% find zero crossings (where F+Psi=0);
q=find(diff(sign(drD)));
q2=q;
for i=1:length(q);
  if (drD(q(i)+1).^2<drD(q(i)).^2)
    q2(i)=q(i)+1;
  end;
end;

subplot(3,2,2*(propnum-1)+1);
plot(gamma2,F,'r',gamma2,-drD,'k',gamma2,-Psi,'b')
hold on;
line([min(gamma2) max(gamma2)],[0 0],'color',[.5 .5 .5]);
set(gca,'xlim',[min(glevels) max(glevels)]);
xlabel('\gamma^n','fontsize',10);
title(['property ',num2str(propnum), ...
	', box ',num2str(boxnum)]);

if (propnum==1)
  plot(gamma2(q2),-drD(q2),'*k')
  disp('Click on zero point for mixing')
  [i,j]=ginput(1);
  if (i>22)
    del=(gamma2(q2)-i).^2;
    del=find(del==min(del));
    q2=q2(del);
  else
    q2=1;
  end;
  drD_ap(1:q2-1,propnum,boxnum)=drD(1:q2-1);
  WA(:,boxnum)=-(F-drD_ap(:,propnum,boxnum));	
	% sign convention: pos.UPWARDS
elseif (propnum==2);
  drD_ap(:,propnum,boxnum)=-[flipud(cumsum(flipud( ...
	diff(WA(:,boxnum).*Spbar) )));0];
else
  drD_ap(:,propnum,boxnum)=-[flipud(cumsum(flipud( ...
        diff(WA(:,boxnum).*Tbar) )));0];
end;

plot(gamma2,-drD_ap(:,propnum,boxnum),'.k');
hold off;
subplot(3,2,2*propnum);
plot(glayers,pfluxrel+pfluxek-diff(F),'b');
hold on;
plot(glayers,pfluxrel+pfluxek ...
	-diff(F-drD_ap(:,propnum,boxnum)),'.k');
hold off;
line([min(gamma2) max(gamma2)],[0 0],'color',[.5 .5 .5]);
set(gca,'xlim',[min(glevels) max(glevels)]);

end;		% loop over properties

if (boxnum<nboxes)
  disp('Paused; hit enter to continue.');
  pause;
  %keyboard;
end;

drD_ap(:,2,boxnum)=0;

dwT=-diff(Tbar).*interp1(gamma2,drD_ap(:,1,boxnum), ...
	glayers,'linear');
drD_ap(:,3,boxnum)=-[flipud(cumsum(flipud(dwT)));0];

end;		% loop over boxes

% define errors on drD_ap
edrD_ap=.5*abs(drD_ap);

eval([' save ' as_dir 'mixing_ap drD_ap edrD_ap; ']);

