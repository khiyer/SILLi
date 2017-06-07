function sill1d_initialize()
% sill1d_initialize
%
% Initialize the environment for the sill model.
%
% Developed by Karthik Iyer, Henrik Svensen and Daniel W. Schmid
%
%
%--------------------------------------------------------------------------

% Check if we need to add paths - look for the gui layout toolbox
pathCell = regexp(path, pathsep, 'split');
if any(strcmpi(fullfile(pwd, '..', 'ext', 'glt', 'layout'), pathCell));
    % Path already set up
    return;
end

% Add paths to required subdirectories
addpath(fullfile(pwd, '..', 'core'));
addpath(fullfile(pwd, '..', 'ext'))
addpath(fullfile(pwd, '..', 'ext', 'glt', 'layout'));
