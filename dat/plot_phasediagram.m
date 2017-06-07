function plot_phasediagram()
% plot_phasediagram
%
% Plot all the phasediagrams contained in the *.mat files in this folder.
%
% Developed by GeoModelling Solutions GmbH
%
% Original author:    Dani Schmid
% Last committed:     $Revision: 0 $
% Last changed by:    $Author: karthik $
% Last changed date:  $Date: 2012-03-07 16:05:55 +0100 (Wed, 07 Mar 2012) $
%--------------------------------------------------------------------------

%% Find all mat files in this folder
listing = dir('*.mat');

%% Plot phase diagrams for H20 and CO2
for i=1:length(listing)
    % Load mat file
    phasediagram    = load(listing(i).name);
    
    % Find name of field that contains the data
    field_name      = fieldnames(phasediagram);
    
    % Plot
    figure('Name', ['Phasediagram: ', listing(i).name(1:end-4)], 'NumberTitle', 'off');
    subplot(121);
    pcolor(phasediagram.(field_name{1}).T, phasediagram.(field_name{1}).P, phasediagram.(field_name{1}).H2O);
    shading interp;
    box on;
    grid on;
    set(gca,'Layer','top');
    xlabel('Temperature [C]');
    ylabel('Pressure [bar]');
    title([listing(i).name(1:end-4), ' Free H_2O [Wt%]']);
    colorbar;
    
    subplot(122);
    pcolor(phasediagram.(field_name{1}).T, phasediagram.(field_name{1}).P, phasediagram.(field_name{1}).CO2);
    shading interp;
    box on;
    grid on;
    set(gca,'Layer','top');
    xlabel('Temperature [C]');
    ylabel('Pressure [bar]');
    title([listing(i).name(1:end-4), ' Free CO_2 [Wt%]']);
    colorbar;
    
    drawnow;
end