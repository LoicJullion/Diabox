function [gvelp,frac,area_tri] = gvelpad(corrpress,gvel,gvel2,press,lat,lon);

% routine to pad gvel2 down to the deepest pressure (from
% corrdepth) for a station pair and then calculate the fraction
% of the area that the bottom velocity estimates represent - 
% this will vary from 1 to zero between the two pressure (from
% corrdepth) estimates.
%modified 08/05/07 HCG: frac needs to be a double array to be written to.

frac = ~isnan(gvel); %needs to be made into a double array
frac=+frac;
gvelp = gvel2;
[pairs] = size(gvel,2);
[area_tri] = ones(size(gvel,2));
interval = press(2) - press(1);
deldistkm      = distance(lat,lon,'km');

for i = 1:pairs;
%establish which is the deepest of the station pair
  if corrpress(i) > corrpress(i+1);
    deep = i;
    shallow = i+1;
  else
    deep = i+1;
    shallow = i;
  end
%find the element numbers for the deepest common level, the
% deepest measurement, and the corrpress of both stations
% for the station pair.
  pdcl = max(find(~isnan(gvel(:,i))));
  pmeas = max(find(~isnan(gvel2(:,i))));
  pdeep = (round((corrpress(deep)-press(1))/interval)) +1;
  pshallow = (round((corrpress(shallow)-press(1))/interval)) +1;
  for j = pmeas+1:pdeep;
    gvelp(j,i) = gvel2(pmeas,i);
  end  
  for j = pdcl+1:pshallow;
   frac(j,i) = 1;
  end
  denominator = pdeep+1-pshallow;
  for j = pshallow+1:pdeep+1;
   frac(j,i) = 1-((j-pshallow)/denominator);
  end
  for j = pdcl+1:pdeep;
   area_tri(i) = area_tri(i)+(frac(j,i)*deldistkm(i)*1000*(sw_dpth(press(j+1),lat(i))-sw_dpth(press(j),lat(i))));
  end
end;

