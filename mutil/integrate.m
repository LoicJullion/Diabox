function yi = integrate_m1(x,y,xi,keyword)

%=========================================================================
% INTEGRATE   4.1   92/11/30   Copyright (C) Phil Morgan 1992
%
% yi = integrate(x,y,xi,{keyword})
%
% DESCRIPTION:
%   Integrates discrete function (x,y) over the intervals xi.
%   NaNs in x or y are automatically stripped but xi in not
%   allowed to contain NaNs.
%
% INPUT:
%   x  = vector of x from (x,y) points. Can have NaNs
%   y  = vector of y from (x,y) points. Can have NaNs
%   xi = vector of intervals to integrate.  No NaNs allowed
%
%   keyword = optional string. Keywords allowed are
%             'view' = view diagnostic plots
%
% OUTPUT:
%   yi = corresponding integration between xi intervals.  
%        thus length(yi) = length(xi)-1.  Also same shape as xi
%        
% AUTHOR:   Phil Morgan.  Copywrite (c) 1992 CSIRO
%
% CALLER:  general purpose
% CALLEE:  none 

%  @(#)integrate.m   Revision: 4.1   92/11/30  Phil Morgan 
%----------------------------------------------------------------
% SET CONSTANTS AND CHECK INPUT ARGUMENTS.  WORK WITH COL VECTORS.
%----------------------------------------------------------------

% modified on 2009.1.20.
% belows are test input file.

%x = [0 2 3 4 6 7];
%y = [10 20 10 30 20 10];
%xi = [1 3 7];

TRUE  = 1;
FALSE = 0;

View  = FALSE;

%keyword = char('view')

if nargin == 4
  View = strcmp(keyword,'view')
end %if

% IF USER PASSED IN ROW VECTOR THEN SET "Transpose" FLAG TO PASS BACK A ROW
[m,n] = size(xi);
%B King at UEA 31 Jul 98. Matlab5 gives awarning if Transpose doesn't exist, so create it.
if(exist('Transpose')~=1) Transpose = [];end;
if m==1 
  xi = xi';
  Transpose = TRUE;   % flag that answer needs to be in row vector
end %if
[mx,nx]=size(x);
if mx == 1
  x = x';
  y = y';
end %if

% USE ONLY VALID X,Y PAIRS.
good = find( ~isnan(x) & ~isnan(y) );
gx   = x( good);
gy   = y( good);
x    = gx;
y    = gy;
   
% SIZE OF X,Y MUST BE THE SAME SO WE GET X,Y PAIRS
if length(x)~=length(y)
  error('INTEGRATE.M: size of x & y without NaNs must be the same')
end %if


% NEED AT LEAST 2 GOOD X,Y PAIRS TO INTERPOLATE "xi" .
if length(x)<2
   disp('INTEGRATE.M: not enought good x,y points need at least 2')
   disp('             Setting integration to NaN')
   yi = NaN*ones(1,length(xi)-1);
   if Transpose
      yi = yi'; 
   end %if 
   return
end %if

% "XI" NOT ALLOWED TO HAVE NANS
if any(any(isnan(xi)))
  error('INTEGRATE.M:  xi has nans!!')
end %if


% XI MUST BE IN THE RANGE OF X
Outrange = FALSE;
toosmall = find(xi < x(1)         );
toobig   = find(xi > x(length(x)) );
if length(toosmall) > 0
   disp('xi smaller than x')
   Outrange = TRUE;
end %if
if length(toobig) > 0
   disp('xi bigger than x')
   Outrange = TRUE;
end %if
if Outrange
    xi
    x1=x(1)
    xn = x(length(x))
   error('INTEGRATE.M: xi must be in range of x ')
end %if

%----------------------------------------------------------------
% BEGIN
%----------------------------------------------------------------
[m,n]=size(x);

% GET MIDPOINT x,y VALUES
delx = x(2:m)-x(1:m-1);
midx = x(1:m-1) + delx/2;
midy = (y(2:m)+y(1:m-1))./2;

% CALC INTEGRAL INTO "tots" AND INTERPOLATE FOR REQUIRED xi. 
tots = cumsum([0; midy.*delx]);
tab  = [x tots];  % tab = 2500x2.

%%%ytot = table1(tab,xi);
% table1 was replaced by below description.
ytot=interp1(tab(:,1),tab(:,2),xi,'linear');

% CALC TOTAL (INTEGRAL)  BEWTEEN THE CONSECUTIVE LIMITS OF XI.
[mm,nn] = size(ytot);
yi      =  ytot(2:mm) - ytot(1:mm-1);
if Transpose
  yi = yi';
end %if

%----------------------------------------------------------------
% DIAGNOSTIC PLOTS IF REQUIRED BY USER
%----------------------------------------------------------------
if View
   hold off
   clg
   subplot(221)
   plot(x,y,'x')
   hold on
   plot(x,y,':')
   hold on
   bar(midx,midy)
   hold on
   plot(xi,zeros(xi),'o')
   disp(' ')
   disp('NOTE: some inaccuracy in bar plots if not evenly spaced.')
   disp('However, results are correct as shown in plot 2')
   disp(' ')
   title('function plotted with integral bewteen points')
   xlabel(' x=data, boxes=integral, o=limits ')
   hold off
   
   plot(xi,ytot,'o')
   grid
   hold on
   plot(xi,ytot)
   title('integrated function')
   hold off
   
   bar(yi)
   title('integration bewteen limits')
   
   disp('hit key to continue...')
   pause
   
end %if
% --------------------------------------------------------------------------