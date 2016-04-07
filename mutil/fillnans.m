function z2=fillnans(z,wrap)
%==========================================================================
% function z2=fillnans(z,wrap)
%
% fills gaps held by NaN's in matrix z via linear
% 2D interpolation to create matrix z2.
% If parameter wrap==1, matrix z is assumed to
% wrap around in the columns [x-dir for z(y,x)].
%
% AUTHOR: Rick Lumpkin
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
qq=find(isfinite(z));
if (length(qq))

    if ~exist('wrap');
      wrap=0;
      disp('Warning: parameter wrap not defined in fillnans.');
    end;

    [n,m]=size(z);
    z2=z;

    if (n>1 & m>1);		% z is a matrix

        % pad edges
        if (wrap)
          z2=[z2(:,m),z2,z2(:,1)];
        else
          z2=[z2(:,1),z2,z2(:,m)];
        end;
        z2=[z2(1,:);z2;z2(n,:)];

        q=find(isnan(z2(2:n+1,2:m+1))==1);
        q=q+2*(floor((q-1)/n)+2)+n-1;

        while (length(q))
          surr=[z2(q+1) z2(q-1) z2(q+n+2) z2(q-n-2) ...
                z2(q+n+3) z2(q+n+1) z2(q-n-1) z2(q-n-3)];
          z2(q)=nanmean(surr')';
          q=find(isnan(z2(2:n+1,2:m+1))==1);
          q=q+2*(floor((q-1)/n)+2)+n-1;
        end;

        z2=z2(2:n+1,2:m+1);

    else	% z is a vector

        needtoflip=0;
        if (n==1)
          z2=z2';
          needtoflip=1;
          n=m;
        end;

        if (wrap)
          z2=[last(z);z2;z(1)];
        else
          z2=[z(1);z2;last(z)];
        end;

        q=find(isnan(z2(2:n+1)))+1;
        while (length(q))
          z2(q)=nanmean([z2(q-1) z2(q+1)]')';
          q=find(isnan(z2(2:n+1)))+1;
        end;

        z2=z2(2:n+1);
        if (needtoflip)
          z2=z2';
        end;


    end

else	% no non-NaN values
    z2=z;
end % First if

return
