function Ln=laynumfield(surf_press,z,bottom);

% function Ln=laynumfield(surf_press,z,bottom);
%
% Calculate a 2D field Ln of dimension
% length(z) by length(bottom), which contains
% layer number for each pixel (as defined by
% surf_press).

Ln=NaN*ones(length(z),length(bottom));

% loop through stations
for j=1:length(bottom);
  % find deepest value of z above bottom
  indmaxz=max(find(z<bottom(j)));
  % loop through z
  for i=1:indmaxz;
    ind=find(surf_press(:,j)>=z(i));
    if length(ind) ind=min(ind);
    else ind=find(surf_press(:,j)== ...
	max(surf_press(:,j)) )+1;
    end;
    if (length(ind))
      Ln(i,j)=ind;
    end;
  end;
end;

