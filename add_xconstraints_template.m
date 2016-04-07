%==========================================================================
%
% add_xconstraints.m
%
%
% Additional property fluxes for sections. (Run "checkgeo" to see
% a map of box and section numbers.)
%
% Constraints actually added by the common m-file addcon.m
%
%  Each call requires the following information:
%   - sectnum: section number
%   - statnum: station pair numbers ([] for all in section)
%   - propnum: the property number (1=mass)
%   - boxnum:  box number
%   - layernum: layer number (nlayers+1 for all)
%   - fluxvalue: constrained value of property flux
%   - fluxerror: error bars on constraint
%
% Because some or all of this information is repeated from one model 
% geometry to the next, with only "boxnum" changing according to the 
% geometry, section-specific m-files with the extension "_con.m" are kept 
% in subdirectory /Path_to_Inversion/inversion_dir/xconstraints/
%
% NOTE: if statnum=[], fluxvalue is assumed to include the Ekman flux.  
% Otherwise, only the geostrophic flux is constrained. This is done because 
% constraints on individual currents (e.g. WBCs) are usually from sub-
% Ekman-layer estimates via current meters, while constraints across an 
% entire section are usually meant to include the Ekman component.
%
% AUTHORS: Rick Lumpkin.
%
% MODIFIED: Loic Jullion.
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
glev;

load dir_loc.mat xconst_dir; % Path to directory constaining the model 
                             % specific constraints
eval(['path(path,''' xconst_dir ''');']); 

%% add constraints on total (ALL boxes) 
% This is to enforce net conservation over all boxes
% property budgets
all_box_cons_switch = 0; % CHANGE 0 TO 1 TO SWITCH ON
if (nboxes>1 & all_box_cons_switch)	% (already done if only one box) 
  for propnum=1:nproperties;
    % assign the A-values for box 1
    AT=A(indb(nlayers+1,propnum,1),:);
    % add the A-values for the other boxes
    for boxnum=2:nboxes;
      AT=AT+A(indb(nlayers+1,propnum,boxnum),:);
    end;
    % add this line to matrix A
    fprintf('%d: ',size(A,1)+1);
    % define bounding sections of entire domain.
    % ASSUMES A SECTION CAN BE SHARED BY TWO AND
    % ONLY TWO BOXES.
    sects=find(sum(geometry,1));
    masserror=sqrt(sum(masssecterror(sects).^2));
    errorratio=masserror./sqrt(sum( ...
	modelerror(nlayers+1,1,:).^2));

    if (propnum==1)
      fprintf('Constraining overall net mass.\n')
      error=masserror;
    elseif (propnum==2)
      fprintf('Constraining overall net salt anom.\n');
      error=errorratio*sqrt(sum( ...
	modelerror(nlayers+1,propnum,:).^2));
    else;
      fprintf('Constraining overall net heat anom.\n');
      error=errorratio*sqrt(sum( ...
        modelerror(nlayers+1,propnum,:).^2));
    end;
    A=[A;zeros(1,size(A,2))];
    A(size(A,1),:)=AT;
    b(length(b)+1)=sum(b(indb(nlayers+1,propnum,:)),3);
    Rnn(length(Rnn)+1,length(Rnn)+1)=error^2;
  end	% loop over properties
end	% if nboxes>1
clear AT sects std_potmp std_salam mean_potmp mean_salam;

% top-to-bottom silica conservation
Sil_cons_switch = 1; % CHANGE 0 TO 1 TO SWITCH ON
if Sil_cons_switch == 1
    for boxnum=1:nboxes;
      if (consSILCAT(boxnum))
        fprintf('%d: Conserving top-to-bottom ',size(A,1)+1);
        fprintf('silica, box %d.\n',boxnum);
        Asil=zeros(1,size(A,2));
        bsil=0;
        sects=find(geometry(boxnum,:));
        for sectnum=1:length(sects);
          [Aline,bvalue]=addsilcat(boxnum,sects(sectnum), ...
        sectlabel,lon,lat,D,A,indE,indb,indv,EkmanE);
          Asil=Asil+Aline;
          bsil=bsil+bvalue;
        end;	% loop over bounding sections of box
        A=[A;Asil]; b=[b;bsil-Si_Input(boxnum)];
        Rnn(length(Rnn)+1,length(Rnn)+1)=(SILCATerror)^2;
      end;	% if consSILCAT(boxnum)
    end; % loop over boxes
    clear Aline Asil Bvalue Bsil;
end

%% model-dependent additional constraints

a05fs_con;    % N flux through Florida Strait
ar16d_con;	% W flux out of Gulf of Cadiz
a08_con;	% N flux across 11S Atl.
