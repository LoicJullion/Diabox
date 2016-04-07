% addcon.m
%
% add the constraints described in addconstraints.m (a model-specific 
% m-file) to A, b and Rnn.
%
% Inputs: sectnum, statnum, propnum, sectnum, boxnum, layernum,
%	fluxvalue, fluxerror.
% 
% Before 22 May 2002:
%   Note that the Ekman flux is included in the constraint
%   only if statnum==[] (i.e. the entire section flux is
%   to be constrained).
% After: the ratio of length(statnum)/(# station pairs in section)
%   is calculated, EkIn is multiplied by this fraction, and 
%   included in the net flux.  The ratio is set to 1 for
%   statnum=[];
%
% AUTHORS: Rick Lumpkin.
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
% The code below won't work for more than one section, if
% a segment of a section is specified by vector statnum.
% Make sure this isn't inadvertantly done by the user.
if (length(statnum) & length(sectnum)>1)
  disp('WARNING: specific station pairs for more than');
  disp('	 one section passed to addcon.m.');
  keyboard
elseif (length(statnum))
  fraction_of_sect=length(statnum)/ ...
	length(find(sectlabel==sectnum));
else
  fraction_of_sect=1;
end;

% calculate vector stns containing indices of section pairs
% matching the section numbers (found in sectnum).
stns=[];
for i=1:length(sectnum);
  q=find(sectlabel==sectnum(i));
  if (length(statnum) & i==1)
    q=q(statnum);
  end;
  stns=[stns;q];
end;

% verify that sectnum and boxnum are consistent
cons=abs(geometry(boxnum,sectnum));
if (length(boxnum)==1)
    % make sure that all sections adjoin box
    if (sum(cons)~=length(cons));
        fprintf('\nWARNING!  Inconsistent constraint:');
        sectnum
        boxnum
        keyboard
    end
else    % multiple boxes, and probably sections too
    if (sum(sum(cons))~=size(cons,1))
        fprintf('\nWARNING! Possible inconsistency?');
        fprintf('see addcon.m')
        sectnum
        boxnum
        keyboard;
    end;
end;

fprintf('%d: ',size(A,1)+1);
if (propnum==1)
  fprintf('Adding constraint: %3.2f +- %3.2f Sv', ...
        fluxvalue/1e6,fluxerror/1e6);
else
  fprintf('Adding non-mass constraint');
end;
fprintf(' across section %d',sectnum(1));
if (length(statnum))
  fprintf(' (pairs %d-%d)',statnum(1),statnum(length(statnum)));
end;
for i=2:length(sectnum);
  fprintf(',%d',sectnum(i));
end;
if (layernum~=nlayers+1)
  fprintf(' in layers %d-%d',min(layernum),max(layernum));
else
  fprintf(' (all layers)');
end;
fprintf(' into box %d.\n',boxnum);
A=[A;zeros(1,size(A,2))];
if (length(layernum)>1 | length(boxnum)>1)
  % multiplicative factors for ref vels
  A(size(A,1),indv(stns))=nansum( ...	% sum over layers
	A(indb(layernum,propnum,boxnum),indv(stns)));
  % multiplicative factors for Ek adjustments
  A(size(A,1),indE(sectnum))=fraction_of_sect*nansum( ...
	A(indb(layernum,propnum,boxnum),indE(sectnum)));
  pfluxrel=nansum(-sum( ...
	D(indb(layernum,propnum,boxnum),indv(stns))')' );
else
  A(size(A,1),indv(stns))= ...
	A(indb(layernum,propnum,boxnum),indv(stns));
  A(size(A,1),indE(sectnum))=fraction_of_sect* ...
	A(indb(layernum,propnum,boxnum),indE(sectnum));
  pfluxrel=-sum( ...
	D(indb(layernum,propnum,boxnum),indv(stns))')' ;
end;

EkmanE(nlayers+1,:,:,:)=nansum(EkmanE(1:nlayers,:,:,:));
EkIn=sum(sum(EkmanE(layernum,propnum,boxnum,sectnum)));
EkmanE(nlayers+1,:,:,:)=[];
% multiply by the fraction of the section specified
% by statnum (fraction_of_sect=1 for statnum=[])
EkIn=EkIn*fraction_of_sect;
if (size(EkIn,4)>1)
    EkIn=sum(EkIn,4);
end;

b(length(b)+1) = fluxvalue-pfluxrel-EkIn;
Rnn(length(Rnn)+1,length(Rnn)+1) = (fluxerror).^2;

clear propnum sectnum boxnum layernum fluxvalue 
clear fluxerror i stns pfluxrel EkIn statnum;

