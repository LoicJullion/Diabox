%==========================================================================
% DIABOX v2.0
% Menu for constructing, solving and analysing box models
%
% DESCRIPTION:
%   Menu for constructing, solving and analysing box model results.
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
% AUTHORS:
%    Loïc Jullion (FSU/MIO). Strongly inspired by Kevin Speer, Rick Lumpkin
%    and the many collaborators who helped me improve my understanding of
%    inverse modelling.
%
%    Rick Lumpkin
%    Kevin Speer
%-----------------
% This is based on the original DOBOX software written by Phil Morgan with 
% many additional revisions by Rick Lumpkin, 2000-2002 and others.
% 
%==========================================================================

box_menu_stop = 0;
box_menu_opt = 1;

while box_menu_stop == 0
    clc
    disp(' ')
    disp('***************************************************************')
    disp('**            DIABOX Version 2.0                             **')
    disp('***************************************************************')
    disp(' 						')
    disp('----- BOX MENU (0) : Model Choices -----')
    disp(' 						')
    disp('      1) Prepare raw Section file for diabox      ')
    disp(' 			                ')   
    disp('      2) Lhs      - Define geometry, build Abase 	')
    disp('      3) Rhs      - Set reflevel, build bbase 	')
    disp('      4) Build    - Define interior diapycnal mixing, build A & b')
    disp(' 				          ')
    disp('      5) Solve inverse model')
    disp('  				')
    disp('      0) Matlab       -    Return to MATLAB 		')
    disp('  							')
    
    defnum = box_menu_opt + 1;
    box_menu_opt = input([' Select a menu number: {' num2str(defnum) '} ']);

  if ~length(box_menu_opt)
    box_menu_opt=defnum;
  end;

  if box_menu_opt == 0   % return to MATLAB
    box_menu_stop = 1;
    
  elseif box_menu_opt == 1   % Process a Section
    disp(' ')
    disp('This gets the raw matlab data file for a section and calculates');
    disp('layer thickness and the average value of the tracer in the areas');
    disp('of the walls of the box and stores the results in a named file');
    disp(' ')
    raw_filename = input('What raw (section) file name ? ','s');

    properties=['mass ';'salam';'potmp']

    sectfile     = dosection(raw_filename,properties);
    clear file
    disp(' ')
    disp('*** Option 1 - dosection - Done ****');
    disp(' ')
    
  elseif box_menu_opt == 2
    file_geom='geo'; % Important to keep the name geo.m for the geometry 
                     % file. Otherwise, the name has to be updated in the
                     % code. Refer to the User guide for more information
                     % on how to construct the geometry of the box.
    eval(file_geom);
    % Make sure the box is well constructed.
    checkgeo('geo');
    disp('Verify the geometry of the box(es):');
    disp('The arrows of the same color should all point inside the same box');
    s = input('Hit "enter" when done','s');

    [Abase, nlayers, nsectpairs] = build_a(geometry, sectfiles, properties);

    [mAbase,nAbase]     = size(Abase);
    [nproperties,m]     = size(properties);
    [nboxes,nsections]  = size(geometry);
    nboxeqns            = (nlayers+1)*nproperties;
    nsurfs              = nlayers + 1;
    clear m
    
    % Define various indices
    [row_beg, row_end, col_beg, col_end, eqn_beg, eqn_end] = ...
          set_index(nboxes,nsections,nproperties,nsectpairs,nlayers,nboxeqns);

    disp(' ')
    disp('*** Option 2 - Build A - Done ****');
    disp(' ')

  elseif box_menu_opt == 3
    % Calculate reference velocities
    % 1) First, define for each section a common value for all the stations 
    % pair 
    refvels(geometry,sectfiles,rv,vmag);
    % 2) Now run "uniquerefvels.m" for sections with particular barotropic 
    %    structure in the initial guess. This file is specific to each
    %    inversion and needs to be modified accordingly.
    uniquerefvels;
    
    % Initialize the matrices containing the a apriori transport
    DL=[]; D=[]; bbase=[]; 
    [DL,D,bbase,reflevel] = rhs_menu(geometry,sectfiles,properties,DL,D,bbase,reflevel);

    disp(' ')
    disp('*** Option 3 - Build DL, D & b -  Done ****');
    disp(' ')

  elseif box_menu_opt == 4
    disp(' ')
    % Get diapycnal velocities
    filename = 'Xcol';
    load dir_loc.mat diapvel_dir;
    [Xcols,XCW] = get_diapvel(filename,nboxes,diapvel_dir);
    
    % We do not had any extra constraints now. They will be added later on.
    % This can probably be removed but it stays here for now.
    Xrows = [];
    
    % CONSTRUCT FULL A AND b MATRIX
    A = [Abase Xcols];
    b = bbase;
	 
    [mA,nA] = size(A);

    % row_beg: First row for each box
    % row_end: Last row for each box
    % col_beg: First col for each ctd section + each Xcol term
    % col_end: Last col for each ctd section + each Xcol term;
    [row_beg,row_end,col_beg,col_end] = ...
          set_xindex(Xrows,Xcols,row_beg,row_end, ...
                     col_beg,col_end,nboxes,nsections);
 
    disp(' ')
    disp('*** Option 4 - set row/col constraints -  Done ****');
    disp(' ')
    disp(' Saving "readyforgauss.mat".');
    load dir_loc.mat output_dir;

    eval(['save ' output_dir 'readyforgauss.mat;']);
  
    % This part is not used for now
%     if (strcmp(input('Recalculate mixing_ap? (y/n) ','s'),'y'))
%       [drD_ap,edrD_ap]=makemixing_ap;
%     else
%       load mixing_ap;
%     end;
%     save readyforgauss;
 
  elseif box_menu_opt == 5
    clear all;
    load dir_loc.mat output_dir;
    matfile_dir = [output_dir 'readyforgauss.mat'];
    if exist(matfile_dir,'file') == 2
       eval(['load ' matfile_dir]);
       gauss;
    else
       warning('The readyforgauss.mat file has not been created yet.');
       warning('Run Step 1 to 4 first.');
    end

    if (~exist('R') | isempty(R))
       disp(' ')
       disp('No solution determined. Returning to main menu')
       disp(' ')
    else

      TF    = zeros(size(A));
      for i=1:mA
        TF(i,:) = A(i,:) .* x';
      end
      clear i
    
      TL = DL;
      TA =TF( row_beg(1):row_end(nboxes), col_beg(1):col_end(nsections))+TL+E;
    end %if

    disp(' ')
    disp('*** Option 6 - method_menu - Done ****');
    disp(' ')   
   
  elseif box_menu_opt == 7
    if (~exist('R') | isempty(R))
       error('dobox: No solution determined in option 6')
    else
       anal_menu;
    end %if
    disp(' ')
    disp('*** Option 7 - anal_menu -  Done ****');
    disp(' ')
   
  else
    box_menu_opt = defnum - 1;

  end %if
end %while

clear box_menu_opt box_menu_stop defnum
return

%--------------------------------------------------------------------


