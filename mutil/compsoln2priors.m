% compsoln2priors
%
% plot the ratio of inversion-derived unknowns 
% (actually their absolute values) to their 
% a priori std.devs, and the noise divided by
% a priori model errors.

% set "apnoise" to 0 to plot n/priors,
% 1 to plot prior noise / prior model errors
apnoise=0;

% create a palette of very light colors
% which will be used to indicate boxes
clf;
colormap default;
if (nboxes>1)
  P=colormap;P=P(1:40,:);
  l1=1:length(P);
  l2=1:(length(P)-1)/(nboxes-1):length(P);
  P2=NaN*ones(nboxes,3);
  for i=1:3;
    P2(:,i)=interp1(l1',P(:,i),l2');
  end;
  P2=1-.25*P2;
  l=rand(size(P2,1),1);
  [l,ind]=sort(l);
  P2=P2(ind,:);
  q=find(max(P2,[],2)-min(P2,[],2)<.2);
  P2(q,:)=.85*P2(q,:);
else
  P2=[.85 .85 .85];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot unknowns divided by a priori %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(211);
maxrat=max(abs(x)./sqrt(diag(Rxx)));
maxrat=max([maxrat 1]);
plot(0,0,'.k','markersize',1);
hold on;

% create shaded regions to indicate boxes
for boxnum=1:nboxes;
  %td=indw(1:nlayers-1,:,boxnum);
  %td=[min(min(td))-.5 max(max(td))+.5];
  %patch([min(td) max(td) max(td) min(td)], ...
	%[0 0 maxrat maxrat],P2(boxnum,:), ...
	%'edgecolor',[.8 .8 .8]);
  % add grid lines to distinguish layers
  for j=min(min(indw(1:nlayers-1,:,boxnum))): ...
        max(max(indw(1:nlayers-1,:,boxnum)));
    plot([j j],[0 maxrat],'color',P2(boxnum,:));
  end;
  %td=indF(1:nlayers-1,:,boxnum);
  %td=[min(min(td))-.5 max(max(td))+.5];
  %patch([min(td) max(td) max(td) min(td)], ...
        %[0 0 maxrat maxrat],P2(boxnum,:), ...
        %'edgecolor',[.8 .8 .8]);
  % add grid lines to distinguish layers
  for j=min(min(indF(1:nlayers-1,:,boxnum))): ...
        max(max(indF(1:nlayers-1,:,boxnum)));
    plot([j j],[0 maxrat],'color',P2(boxnum,:));
  end;

end;

% plot vref
plot(abs(x(indv))./sqrt(diag(Rxx(indv,indv))), ...
        'color',[.5 .5 .5]);
tiks=(1:length(indv))';
% loop through boxes, plotting w*s for each property
for boxnum=1:nboxes;
  for propnum=1:nproperties;
    td=indw(1:nlayers-1,propnum,boxnum);
    plot(td,abs(x(td))./sqrt(diag(Rxx(td,td))), ...
        'color',propcolr(propnum,:));
    tiks=[tiks;(1:length(td))'];
  end;
  td=indF(1:nlayers-1,1,boxnum);
  plot(td,abs(x(td))./sqrt(diag(Rxx(td,td))), ...
        '.','color',[1 .2 .2], ...
        'markersize',6);
  tiks=[tiks;(1:length(td))'];
  td=indF(1:nlayers-1,2,boxnum);
  plot(td,abs(x(td))./sqrt(diag(Rxx(td,td))), ...
        '.','color',[.5 .5 1], ...
        'markersize',6);
  tiks=[tiks;(1:length(td))'];
end;
% plot E* for each section
td=indE;
plot(td,abs(x(td))./sqrt(diag(Rxx(td,td))),'*k');
tiks=[tiks;(1:length(td))'];

if (maxrat>1);
  plot([1 length(x)],[1 1],':k');
end;
hold off;
title('x / priors');
axis([1,length(x),0,maxrat]);
tk=(1:length(x))';
q=find((tiks==5*round(tiks/5) & tk>length(indv)) | ...
	tk>=min(indE));
set(gca,'xtick',tk(q), ...
  'xticklabel',num2str(tiks(q)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot a posteriori noise for each conservation %
% equation vs. the a proiri model error         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (apnoise)
  nplot=b;
else
  nplot=n;
end;
q=find(nplot==0);
nplot(q)=NaN;

subplot(212);
plot(0,0,'.k','markersize',1);
hold on;
maxrat=max(abs(nplot)./sqrt(diag(Rnn)));
maxrat=max([maxrat 1]);
minrat=min(abs(nplot)./sqrt(diag(Rnn)));

% create shaded regions to indicate boxes
for boxnum=1:nboxes;
  %td=indb(1:nlayers,:,boxnum);
  %td=[min(min(td))-.5 max(max(td))+.5];
  %patch([min(td) max(td) max(td) min(td)], ...
        %[0 0 maxrat maxrat],P2(boxnum,:), ...
        %'edgecolor',[.8 .8 .8]);
  % add grid lines to distinguish layers
  for j=min(min(indb(1:nlayers,:,boxnum))): ...
        max(max(indb(1:nlayers,:,boxnum)));
    plot([j j],[1e-10 maxrat],'color',P2(boxnum,:));
  end;
end;

% plot n/priors
tiks=[];
for boxnum=1:nboxes;
  for propnum=1:nproperties;
    % layer-by-layer model errors/priors
    td=indb(1:nlayers,propnum,boxnum);
    plot(td,abs(nplot(td))./sqrt(diag(Rnn(td,td))), ...
        'color',propcolr(propnum,:));
    tiks=[tiks;(1:length(td))'];
    % vertically-integrated model errors/prior
    td=indb(nlayers+1,propnum,boxnum);
    plot(td-nlayers/2,abs(nplot(td))./sqrt(diag(Rnn(td,td))), ...
        '.','color',propcolr(propnum,:), ...
        'markersize',20);
    tiks=[tiks;(1:length(td))'];
  end;
end;
% errors in meeting added constraints/priors
td=max(max(max(indb)))+1:length(nplot);
if (length(td))
  plot(td,abs(nplot(td))./sqrt(diag(Rnn(td,td))), ...
    'ok','markersize',6);
  tiks=[tiks;td'];
end;
if (maxrat>1);
  plot([1 length(nplot)],[1 1],':k');
end;
hold off
if (apnoise)
  title('a priori noise / model errors')
else
  title('n / priors');
end;
if (apnoise)
  axis([1,length(nplot),minrat,maxrat]);
  %set(gca,'yscale','log')
else
  axis([1,length(nplot),0,maxrat]);
end;
tk=(1:length(nplot))';
q=find(tiks==5*round(tiks/5) |  ...
	tk>max(max(max(indb))));
set(gca,'xtick',tk(q), ...
  'xticklabel',num2str(tiks(q)), ...
  'tickdir','out','ticklength',[.005 .001]);

clear nplot

