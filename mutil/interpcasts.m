function [A,LayP]=interpcasts(sp1,sp2,bottom,dx,nsurfs)
%==========================================================================
% function [A,LayP]=interpcasts(sp1,sp2,bottom,dx,nsurfs);
%
% Calculate A(nsurfs), a vector of areas (units of [bottom]*[dx]) for each 
% layer as defined by level pressure values contained in sp1 and sp2 (the 
% latter located distance dx from sp1 cast).
% Also calculate the mean pressure of each level.
%
% AUTHORS: Rick Lumpkin, 29 May 2002
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
x=0:.01:1;
z=0:.01:1;
Lay=NaN*ones(length(z),length(x));

za1=[sp1(find(isfinite(sp1)))];
za2=[sp2(find(isfinite(sp2)))];
xa1=[zeros(size(sp1(find(isfinite(sp1)))))];
xa2=[ones(size(sp2(find(isfinite(sp2)))))];
za1=[sp1(find(isfinite(sp1)))];
za2=[sp2(find(isfinite(sp2)))];
maxz=max([max(za1) max(za2)]);
za1=za1/maxz;za2=za2/maxz;
botz=[bottom(1) bottom(2)]/maxz;
Laya1=[find(isfinite(sp1))];
Laya2=[find(isfinite(sp2))];
botz=interp1([0 1],botz,x,'linear');

Lay(:,1)=fillnans(interp1(za1,Laya1,z)',0);
Lay(:,size(Lay,2))=fillnans( ...
	interp1(za2,Laya2,z)',0);

for jj=1:length(z);
  Lay(jj,2:length(x)-1)=interp1( ...
	[0 1],[Lay(jj,1) Lay(jj,length(x))], ...
	x(2:length(x)-1),'linear');
end;
for jj=1:length(x);
  Lay(find(z>botz(jj)),jj)=NaN;
end;

z=z*max(bottom);
x=x*dx;

% define mean level pressures
zmat=z'*ones(1,length(x));
LayP=NaN*ones(nsurfs+2,1);
for i=1:length(LayP);
  q=find(Lay>i-.5 & Lay<i+.5);
  if (length(q))
    LayP(i)=mean(zmat(q));
  end;
end;
LayP(find(LayP==min(LayP)))=0;
LayP(find(LayP==max(LayP)))=max(bottom);

% interpolate missing interior values
% (can happen with extremely thin layers)
q=find(isfinite(LayP));
q=min(q):max(q);
lp=LayP(q);
q1=find(isnan(lp));
q2=find(isfinite(lp));
if (length(q1))
  lp(q1)=interp1(q2,lp(q2),q1,'linear');
  LayP(q)=lp;
end;

% convert to layer number
Lay=floor(Lay);
%patchmap2(x/1e3,-z,Lay);colorbar;
%hold on;H=contour(x/1e3,-z,Lay,1:nsurfs+1);
%clabel(H);hold off;

% calculate total area, including
% bottom triangle
qtot=find(isfinite(Lay));
Atot=length(qtot)*abs(median(diff(x)) ...
	*median(diff(z)));

zmat=z'*ones(1,length(x));

A=NaN*ones(nsurfs+1,1);
for i=1:length(A);
  q=find(Lay==i);
  A(i)=Atot*length(q)/length(qtot);
end;
return