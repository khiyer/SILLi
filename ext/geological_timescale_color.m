function [COLOR_PERIOD, COLOR_STAGE] = geological_timescale_color(Age, color_period, color_stage)
% geological_timescale_color
% 
% Returns the period color for given ages.
%
% Original author:    Schmid
% Last committed:     $Revision: 243 $
% Last changed by:    $Author: schmid $
% Last changed date:  $Date: 2012-03-07 16:05:55 +0100 (Wed, 07 Mar 2012) $
%--------------------------------------------------------------------------

% Nargin
if nargin<1
    errordlg('Not enough input');
end
if nargin<2
    color_period    = 1;
end
if nargin<3
    color_stage     = 1;
end

% Get Data
[Period_age, Period_name, Period_rgb, Stage_age, Stage_name, Stage_rgb] = geological_timescale_data();

% Period Based Colors
if color_period
    % Initialize Output
    COLOR_PERIOD   = zeros(length(Age), 3);
    
    % Assign Colors
    for k=1:length(Age)
        for l=1:length(Period_name)
            if (Period_age(l+1)>=Age(k)) && (Age(k)>=Period_age(l))
                COLOR_PERIOD(k, :)   = Period_rgb(l, :);
                break;
            end
        end
    end
else
    COLOR_PERIOD   = [];
end


% Stage Based Colors
if color_stage
    % Initialize Output
    COLOR_STAGE   = zeros(length(Age), 3);
    
    % Assign Colors
    for k=1:length(Age)
        for l=1:length(Stage_name)
            if (Stage_age(l+1)>=Age(k)) && (Age(k)>=Stage_age(l))
                COLOR_STAGE(k, :)   = Stage_rgb(l, :);
                break;
            end
        end
    end
else
    COLOR_STAGE   = [];
end