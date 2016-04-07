% This is a template geometry file which is copied inside every new
% inversion directory creating using the mfile "create_new_inversion.m"
%
% This is the file where many important features of the inverse model are
% set. It needs to be build with care.
%
% Loic - FSU - January 2014
%--------------------------------------------------------------------------
%% Define the sections you want to use for the inversion

% 1) Chose the directory
load dir_loc.mat sect_dir; % Load the path to the directory containing the 
                           % sections. dir_loc.mat was created by the 
                           % create_new_inversion.m script
                           
% 2) Enter the name of all the sections you want to use.
% In this example we have 2 boxes in the tropical and north atlantic
sectfiles=[	sect_dir,'a08.mat    '; ...     
            sect_dir,'a05_fs.mat '; ...
            sect_dir,'a05int.mat ';...
            sect_dir,'ar16d.mat  ';...
            sect_dir,'ar19a.mat  '];
%% Now define the geometry. 
% See the documentation on how to create the correct geometry file.
% Here is a simple example with 2 boxes in the Atlantic
% Initialize your geometry file
geometry=zeros(1,size(sectfiles,1)); bn=1; %bn = box number
% Fill in 1s for the sections corresponding to the first box.
% CAtl 11S to 24N
geometry(bn,sectname2sectnum('a05int',sectfiles))=-1;
geometry(bn,sectname2sectnum('a05fs',sectfiles))=-1;
geometry(bn,sectname2sectnum('a08',sectfiles))=1;
% NAtl 24N to 48N
bn=bn+1;
geometry(bn,sectname2sectnum('ar19a',sectfiles))=-1;
geometry(bn,sectname2sectnum('ar16d',sectfiles))=1;
geometry(bn,sectname2sectnum('a05int',sectfiles))=1;
geometry(bn,sectname2sectnum('a05fs',sectfiles))=1;

% If you want more boxes, simply copy and paste after bn = bn+1 and add new
% raws to the variable geometry
%% Set some default values for the reference velocities, level and errors: 0
% reference velocities referenced at the bottom (-1) and 1 cm/s error.
% All these default values can be changed further along in DIABOX (see
% documentation)

% Define some default settings for the reference velocities, reference
% levels and uncertainties on reference velocities
rv       = zeros(1,size(geometry,2));
reflevel = -1*ones(size(rv));
vmag     = .01*ones(size(rv));

%% Properties to include in the inversion. 
% The minimum is mass, sal and temp but other tracers can be added.
properties=['mass '; 'salam';'potmp'];

%% Define default uncertainties on the transport
[masssecterror,masstotalerror]= definemasserrors(sectfiles,geometry);

%% Silicate conservation in the different boxes.
% consSILCAT(nboxes): set to 1 to conserve net silicate 
% Default is set to 0
consSILCAT=zeros(1,size(geometry,1));
%consSILCAT(1:10)=1;
