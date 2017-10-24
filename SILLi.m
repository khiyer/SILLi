% SILLi
%
% Main run file for 1D sill model.
%
% Developed by Karthik Iyer, Henrik Svensen and Daniel W. Schmid
%
%--------------------------------------------------------------------------
clear all
codeloc                 = fileparts(mfilename('fullpath')); % Command only works within script
% codeloc                 = pwd; % Command does not work if script is added
% to path

%% User Input
% Input file
Filename                = fullfile(codeloc, 'examples', '1d_sill_input_dyke.xlsx');

% Resolution - specified through minimum grid spacing and minimum number of
% points inside sedimentary and sill layers. Whichever measure produces
% higher resolution is taken.
resolution.dz_sed       = 5;            % minimum spacing in sediments (m)
resolution.dz_sill      = 1;            % minimum spacing in sills (m)
resolution.pts_sed      = 20;           % minimum number of points in sediments
resolution.pts_sill     = 50;           % minimum number of points in sills

%% Initialize (Paths)
disp('SILLi Model');
disp(' Initializing');
run(fullfile(codeloc, 'core', 'sill1d_initialize'));

%% Check If Already Run
% If there exists a mat file with the same filename as the Excel input file
% then the possibility is offered to load it and bypass the computation

% Check if mat file with same filename as input exists
Ind         =  strfind(Filename, '.');
Filename_mat= [Filename(1:Ind(end)), 'mat'];
if exist(Filename_mat, 'file')
    disp(' Found result file');
    Choice  = questdlg({'Result file for this input found.', 'If you did not change the input file or the resolution you can bypass running the model again.'}, 'Load or Run?', 'Load', 'Run', 'Load');
    if isempty(Choice)
        error('sill1d_main:no_choice', 'You must decide to load previous results or run the model.');
    end
else
    Choice  = 'Run';
end

%% Run Model
if strcmp(Choice, 'Run') 
    %% - Read Input File
    disp(' Reading input file');
    [rock, sill, fluid, welldata]   = sill1d_input(Filename);
    
    %% - Show input
    sill1d_gui_handle = sill1d_analyze(rock, sill, welldata);
    
    %% - Compute
    disp(' Running numerical model');
    [result, release] = sill1d_compute(rock, sill, fluid, welldata, resolution, codeloc);
    
    %% - Save
    disp(' Saving model results');
    save(Filename_mat, 'rock', 'sill', 'welldata', 'result', 'release');
    
    %% - Close Input Visualization
    delete(sill1d_gui_handle);
else
    disp(' Loading model results');
    load(Filename_mat);
end

%% Analyze Results
disp(' Showing model results');
sill1d_analyze(rock, sill, welldata, result, release, codeloc);
