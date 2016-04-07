
%... This mfile sets up the interfacial statements
%... for the interior diapycnal fluxes in DOBOX

% input the geometry file
disp('Loading geometry file "geo".')
geo; 
nboxes=size(geometry,1);
clear directory geometry properties sectfiles;

% directory containing mprop_atl* files
dir='';

%... load mprop_atl*.mat files
disp('Hard-coded: conserving mass, salam, temp')

for i = 1:nboxes 
  command =['mprop_atl',int2str(i), ' area* msalt* mtemp*'];

  % load data
  eval(['load ' dir command]);
  % define salinity anomaly
  eval(['msaltan_atl',int2str(i), ...
	' = (msalt_atl',int2str(i), ' - 35)*1e-3;']);

%... use diag to form the diagonal matrix with areas on the 
%... diagonal and -area on the next lower diagonal.
  
  eval(['m = length(area_atl', int2str(i), ');']);
  eval(['[lg,icolg] = lastgood(mtemp_atl', int2str(i), ',''' '0' ''' );'])

  eval(['mcol_atl', int2str(i), ' = diag(area_atl', int2str(i), ') +' ...
        'diag(-area_atl', int2str(i), '(1:m-1), -1);'])
  
  eval(['tcol_atl', int2str(i), ' = diag(area_atl', int2str(i), ...
   '.*mtemp_atl', int2str(i), ') + diag(-area_atl', int2str(i), ...
        '(1:m-1) .*mtemp_atl', int2str(i), '(1:m-1), -1);'])
  
  eval(['scol_atl', int2str(i), ' = diag(area_atl', int2str(i), ...
   '.*msaltan_atl', int2str(i), ') + diag(-area_atl', int2str(i), ...
        '(1:m-1) .*msaltan_atl', int2str(i), '(1:m-1), -1);'])
  
%... add two rows of zeros to end of each diagonal matix. These 
%... balance the xcol to give 18 rows for each property flux.
%... The last good value of the xcol is the flux across the bottom
%... boundary. For the Pacific this is at 895 dbars.

  eval(['mcol_atl', int2str(i), ' = ' ...
     '[mcol_atl', int2str(i), '; zeros(2,m)];']);
   
  eval(['tcol_atl', int2str(i), ' = ' ...
     '[tcol_atl', int2str(i), '; zeros(2,m)];']);   
  
  eval(['scol_atl', int2str(i), ' = ' ...
     '[scol_atl', int2str(i), '; zeros(2,m)];']);

%... Form xcol matrix for each box

  eval(['xcol_atl0 = zeros(size(mcol_atl', ...
	num2str(i),'));']);

  eval(['Xcols',int2str(i),' = ', ...
	'[mcol_atl',num2str(i),' xcol_atl0 xcol_atl0;', ...
	' xcol_atl0 scol_atl',num2str(i),' xcol_atl0;', ...
	' xcol_atl0 xcol_atl0 tcol_atl',num2str(i),'];']);


end;		% loop over boxes

% Form overall Xcols matrix 
Xcols=Xcols1;
n=size(Xcols1,1);m=size(Xcols1,2);
for i=2:nboxes
  eval(['Xcols=[Xcols,zeros((i-1)*n,m);', ...
	'zeros(n,(i-1)*m),Xcols',num2str(i),'];']);
end;

Xcol_int = change(Xcols,'==',NaN,0);

[nr,nc] = size(Xcols);

XCW_int =  ones(1,nc);

Xcols=Xcol_int;
XCW=XCW_int;

save Xcol Xcols XCW

