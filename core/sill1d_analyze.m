function sill1d_gui_handle = sill1d_analyze(rock, sill, welldata, result, release, codeloc)
% sill1d_analyze
%
% Postprocessing of sill model results.
%
% Developed by Karthik Iyer, Henrik Svensen and Daniel W. Schmid
%

%--------------------------------------------------------------------------

yr              = 365*24*60*60;

%% Initialize Figure
Screensize 	= get(0, 'ScreenSize');
x_res     	= Screensize(3);
y_res     	= Screensize(4);

fracx     	= .6;
fracy      	= .6;
gui_width  	= round(x_res*fracx);
gui_height	= round(y_res*fracy);
gui_x      	= round((x_res-gui_width)/2);
gui_y      	= round((y_res-gui_height)/2);

sill1d_gui_handle = figure(...
    'Units', 'pixels','pos',round([gui_x gui_y gui_width  gui_height]),...
    'Name', 'SILLi',...
    'Tag', 'sill1d_gui_handle',...
    'NumberTitle', 'off', ...
    'DockControls', 'off', ...
    'Toolbar', 'none', ...
    'Units', 'Pixels');

toolbar = uitoolbar('parent', sill1d_gui_handle, 'handleVisibility', 'off', 'tag', 'FigureToolBar');

SaveFigure  = uitoolfactory(toolbar, 'Standard.SaveFigure');
set(SaveFigure, 'ClickedCallback', @(a,b) save_hires_fig, 'tooltipstring', 'Save Figure')
ZoomIn      = uitoolfactory(toolbar, 'Exploration.ZoomIn');
set(ZoomIn, 'Separator', 'on');
uitoolfactory(toolbar, 'Exploration.ZoomOut');
uitoolfactory(toolbar, 'Exploration.Pan');
DataCursor  = uitoolfactory(toolbar, 'Exploration.DataCursor');
set(DataCursor, 'Separator', 'on');
uitoolfactory(toolbar, 'Exploration.Brushing');

%% Main Layout - 3 Tabs
tab_panel           = uix.TabPanel('Parent', sill1d_gui_handle, 'Padding', 5);
h_tp_input      	= uix.HBox('Parent', tab_panel);
h_tp_result         = uix.VBox('Parent', tab_panel);
h_tp_release    	= uix.HBox('Parent', tab_panel);
tab_panel.TabTitles = {'Input', 'Results', 'Release'};
tab_panel.Selection = 2;

%% Input Visualization
%% - Layout
h_uc            = uicontainer('parent', h_tp_input);
ax1_1           = axes('parent', h_uc, 'ActivePositionProperty', 'outerposition');
h_uc            = uicontainer('parent', h_tp_input);
ax1_2           = axes('parent', h_uc, 'ActivePositionProperty', 'outerposition');
h_uc            = uicontainer('parent', h_tp_input);
ax1_3           = axes('parent', h_uc, 'ActivePositionProperty', 'outerposition');
h_uc            = uicontainer('parent', h_tp_input);
ax1_4           = axes('parent', h_uc, 'ActivePositionProperty', 'outerposition');
h_uc            = uicontainer('parent', h_tp_input);
ax1_5           = axes('parent', h_uc, 'ActivePositionProperty', 'outerposition');

%% -- Link Axes
linkaxes([ax1_1, ax1_2, ax1_3, ax1_4, ax1_5], 'y');

%% - Lithology
min_ax = -min([rock.Tops; sill.Tops])./1e3;
max_ax = -max([rock.Tops; sill.Tops+sill.Thick])./1e3;
axis(ax1_1, [0 1 floor(max_ax) ceil(min_ax)]);
axes(ax1_1); % Needed because hatchfill does not have parent specification
hold(ax1_1, 'on');
% Plot sills and sedimentary rocks as patches
for i = 1:length(rock.Tops)-1
    % Geological Color
    [~,A] = geological_timescale_color(rock.Ages(i));
    
    % Patches
    if i<length(rock.Tops)-1
        h_p     = patch([0 1 1 0], -[rock.Tops(i+1) rock.Tops(i+1) rock.Tops(i) rock.Tops(i)]./1e3, A, 'Parent',ax1_1);
    else
        h_p     = patch([0 1 1 0], -[max([rock.Tops; sill.Tops+sill.Thick]) max([rock.Tops; sill.Tops+sill.Thick]) rock.Tops(i) rock.Tops(i)]./1e3, A, 'Parent',ax1_1);
    end
    
    % Add contextmenu to patch that shows age
    uicm    = uicontextmenu;
    h_p.UIContextMenu = uicm;
    if isnan(rock.Ero_t(i))
        uimenu(uicm, 'Label', [rock.Name{i} ' [' num2str(rock.Ages(i)) ' Ma]']);
    else
        uimenu(uicm, 'Label', [rock.Name{i} ' [' num2str(rock.Ages(i)) ' Ma | ' num2str(rock.Ero_t(i)) ' Ma]']);
        % Add erosion pattern over areas that are eroded
        hatchfill(h_p, 'cross', 90, 20, A);
    end
end

% Add speckle over the sills
for i = 1:length(sill.Tops)
    p = patch( [0 1 1 0], -[sill.Tops(i)+sill.Thick(i) sill.Tops(i)+sill.Thick(i) sill.Tops(i) sill.Tops(i)]./1e3, 'k', 'Parent',ax1_1);
    hatchfill(p, 'speckle');
end

% Ornaments
set(ax1_1, 'XTick', [], 'XTickLabel', [], 'Box', 'on');
ylabel(ax1_1, 'TVDSS [km]');
title(ax1_1, 'Lithological Column');

%% - Density
h_p = plot(rock.Rho, -rock.Tops./1e3, 'o-', 'Color', [0 102/255 204/255], 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Parent', ax1_2); % Blue
add_menu(h_p);
ylabel(ax1_2, 'TVDSS [km]')
title( ax1_2, 'Density [kg/m^3]')
set(   ax1_2, 'xaxisLocation', 'top')
grid(  ax1_2, 'on')
ylim(  ax1_2, [floor(max_ax) ceil(min_ax)]);

%% - Porosity
h_p = plot(rock.Phi, -rock.Tops./1e3, 'o-', 'Color', [255/255 128/255 0], 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Parent',ax1_3); % Orange
add_menu(h_p);
ylabel(ax1_3, 'TVDSS [km]')
title( ax1_3, 'Porosity [fraction]')
set(   ax1_3, 'xaxisLocation', 'top')
grid(  ax1_3, 'on')
ylim(  ax1_3, [floor(max_ax) ceil(min_ax)]);

%% - TOC
h_p = plot(rock.Toc.*100, -rock.Tops./1e3, 'o-', 'Color', [0 204/255 102/255], 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Parent',ax1_4); % Green
add_menu(h_p);
ylabel(ax1_4, 'TVDSS [km]')
title( ax1_4, 'TOC Content [wt%]')
set(   ax1_4, 'xaxisLocation', 'top')
grid(  ax1_4, 'on')
ylim(  ax1_4, [floor(max_ax) ceil(min_ax)]);

%% - Conductivity
h_p = plot(rock.K, -rock.Tops./1e3, 'o-', 'Color', [255/255 51/255 51/255], 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Parent',ax1_5); % Red
add_menu(h_p);
ylabel(ax1_5, 'TVDSS [km]')
title( ax1_5, 'Conductivity [W m^-^1 K^-^1]')
set(   ax1_5, 'xaxisLocation', 'top')
grid(  ax1_5, 'on')
ylim(  ax1_5, [floor(max_ax) ceil(min_ax)]);

%% Results Visualization
%% - Check if we have data
if ~exist('result', 'var') || isempty(result)
    tab_panel.Selection = 1;
    drawnow;
    return;
end

%% - Number of timesteps
no_tstep    = length(result.Time);

%% - Layout
h_t2_hb1	= uix.HBox('Parent', h_tp_result);
h_t2_hb2	= uix.HBox('Parent', h_tp_result);
h_tp_result.Heights = [20, -1];

%% -- Play buttons
PlayButtons = load(fullfile(codeloc, 'ext', 'PlayButtons.mat'));

%   START button
obj.sill1d_play_go_start = ...
    uicontrol('Parent', h_t2_hb1, 'style', 'pushbutton',...
    'cdata', PlayButtons.Start, 'units', 'pixels',...
    'tag', 'sill1d_play_go_start',...
    'callback',  @(a,b)  sill1d_play);

%   PLAY BACK button
obj.sill1d_play_backward = ...
    uicontrol('Parent', h_t2_hb1, 'style', 'pushbutton',...
    'cdata', PlayButtons.PlayBack, 'units', 'pixels',...
    'tag', 'sill1d_play_backward',...
    'callback',  @(a,b)  sill1d_play);

%   STOP button
obj.sill1d_play_stop = ...
    uicontrol('Parent', h_t2_hb1, 'style', 'pushbutton',...
    'cdata', PlayButtons.Stop, 'units', 'pixels',...
    'tag', 'sill1d_play_stop',...
    'callback',  @(a,b)  sill1d_play);

%   PLAY button
obj.sill1d_play_forward = ...
    uicontrol('Parent', h_t2_hb1, 'style', 'pushbutton',...
    'cdata', PlayButtons.Play, 'units', 'pixels',...
    'tag', 'sill1d_play_forward',...
    'callback',  @(a,b)  sill1d_play);

%   END button
obj.sill1d_play_go_end = ...
    uicontrol('Parent', h_t2_hb1, 'style', 'pushbutton',...
    'cdata', PlayButtons.End, 'units', 'pixels',...
    'tag', 'sill1d_play_go_end',...
    'callback',  @(a,b)  sill1d_play);

%   Time
obj.sill1d_time = ...
    uicontrol('Parent', h_t2_hb1, 'style', 'edit', 'String', '65.5', ...
    'tag', 'sill1d_time', ...
    'callback',  @(a,b)  sill1d_play);
obj.sill1d_time_info = ...
    uicontrol('Parent', h_t2_hb1, 'style', 'text', 'String', '[Ma] (1/100)', ...
    'HorizontalAlignment', 'left', ...
    'tag', 'sill1d_time_info');

% Empty for the rest
uix.Empty('Parent', h_t2_hb1);

% Sizes
h_t2_hb1.Widths     = [repmat(20, 1, 5), 75, 100, -1];

%% -- Axes
h_uc            = uicontainer('parent', h_t2_hb2);
ax.res(1)       = axes('parent', h_uc, 'ActivePositionProperty', 'outerposition');
h_uc            = uicontainer('parent', h_t2_hb2);
ax.res(2)       = axes('parent', h_uc, 'ActivePositionProperty', 'outerposition');
h_uc            = uicontainer('parent', h_t2_hb2);
ax.res(3)       = axes('parent', h_uc, 'ActivePositionProperty', 'outerposition');
h_uc            = uicontainer('parent', h_t2_hb2);
ax.res(4)       = axes('parent', h_uc, 'ActivePositionProperty', 'outerposition');
h_uc            = uicontainer('parent', h_t2_hb2);
ax.res(5)       = axes('parent', h_uc, 'ActivePositionProperty', 'outerposition');
h_uc            = uicontainer('parent', h_t2_hb2);
ax.res(6)       = axes('parent', h_uc, 'ActivePositionProperty', 'outerposition');
h_uc            = uicontainer('parent', h_t2_hb2);
ax.res(7)       = axes('parent', h_uc, 'ActivePositionProperty', 'outerposition');

%% -- Link Axes
linkaxes(ax.res, 'y');

%% - Plot Present Day
plot_results(no_tstep);

%% Release Visualization
%% - Layout
h_uc            = uicontainer('parent', h_tp_release);
ax3_1           = axes('parent', h_uc, 'ActivePositionProperty', 'outerposition');
h_uc            = uicontainer('parent', h_tp_release);
ax3_2           = axes('parent', h_uc, 'ActivePositionProperty', 'outerposition');

%% - CO2 Release
% Organic
x2              = [];
Emp_times       = unique(sill.E_time);
for jj = 1:length(Emp_times)
    if jj == 1
        sill_t  = max(sill.Thick(sill.E_time==Emp_times(jj)))/2;
        s_time  = ceil(10*sill_t^2*max(rock.Rho)*max(rock.Cp)/min(rock.K)/yr);
        x2      = unique([0:5000:Emp_times(jj)*1e6-s_time Emp_times(jj)*1e6-s_time:1:Emp_times(jj)*1e6+100]);
    else
        x2      = unique([x2 unique([x2(end):5000:Emp_times(jj)*1e6-s_time Emp_times(jj)*1e6-s_time:1:Emp_times(jj)*1e6+100])]);
    end
end
x2              = unique([x2 x2(end):5000:round(max(-(result.Time-max(result.Time)))) round(max(-(result.Time-max(result.Time))))]);                
y2              = interp1(-(result.Time-max(result.Time))./1e6, sum(release.CO2_org), x2./1e6);     

[h_ax, h_p1, h_p2] = plotyy(...
    -(result.Time-max(result.Time))./1e6, sum(release.CO2_org)./1e3, ...
    x2./1e6, [0 abs(diff(y2))./diff(x2)], 'Parent', ax3_1);
h_ax(2).YLim(1) = 0;
h_ax(1).XLim(2) = round(max(-(result.Time-max(result.Time))./1e6));
h_ax(2).XLim(2) = round(max(-(result.Time-max(result.Time))./1e6));

% Which axis is on top
ax_ontop = 2; %
uistack(h_ax(ax_ontop));
set(h_ax(ax_ontop), 'Color', 'none', 'Hittest', 'off');
set(h_ax(setdiff(1:2, ax_ontop)), 'Color', 'w');
%-------------------

add_menu(h_p1);
add_menu(h_p2);
ylabel(h_ax(1), 'Organic CO_2 Released [ton/m^2]');
ylabel(h_ax(2), 'Organic CO_2 Released [kg/m^2/yr]');
title( ax3_1, 'Organic CO_2 Released');
xlabel(ax3_1, 'Time [Ma]');
set(  h_ax(1), 'xdir', 'reverse');
set(  h_ax(2), 'xdir', 'reverse');
grid(  h_ax(1), 'on');

%% - Inorganic
y3              = interp1(-(result.Time-max(result.Time))./1e6, sum(release.CO2_rel), x2./1e6);

[h_ax, h_p1, h_p2] = plotyy(...
    -(result.Time-max(result.Time))./1e6, sum(release.CO2_rel)./1e3, ...
     x2./1e6, [0 abs(diff(y3))./diff(x2)], 'Parent', ax3_2);
h_ax(2).YLim(1) = 0;
h_ax(1).XLim(2) = round(max(-(result.Time-max(result.Time))./1e6));
h_ax(2).XLim(2) = round(max(-(result.Time-max(result.Time))./1e6));

% Which axis is on top
ax_ontop = 2; %
uistack(h_ax(ax_ontop));
set(h_ax(ax_ontop), 'Color', 'none', 'Hittest', 'off');
set(h_ax(setdiff(1:2, ax_ontop)), 'Color', 'w');
%-------------------

add_menu(h_p1);
add_menu(h_p2);
ylabel(h_ax(1), 'Inorganic CO_2 Released [ton/m^2]');
ylabel(h_ax(2), 'Inorganic CO_2 Released [kg/m^2/yr]');
title( ax3_2, 'Inorganic CO_2 Released');
xlabel(ax3_2, 'Time [Ma]');
set(  h_ax(1), 'xdir', 'reverse');
set(  h_ax(2), 'xdir', 'reverse');
grid(  h_ax(1), 'on');

%% fun sill1d_play
    function sill1d_play
        % Who is calling
        wcbo            = gcbo;
        Whoiscalling    = get(wcbo, 'tag');
        
        % Get currently plotted tstep
        tstep   = getappdata(sill1d_gui_handle, 'tstep');
        
        switch Whoiscalling
            case 'sill1d_play_go_start'
                plot_results(max(2, tstep-1));
                
            case 'sill1d_play_backward'
                
                % Set an interrupt flag and write into storage
                % flag = 0;
                % setappdata(sill1d_gui_handle, 'flag', flag);
                
                %  Enable appropreate play buttons
                set(obj.sill1d_play_go_start,   'enable', 'off');
                set(obj.sill1d_play_backward,   'enable', 'off');
                set(obj.sill1d_play_stop,       'enable', 'on');
                set(obj.sill1d_play_forward,    'enable', 'off');
                set(obj.sill1d_play_go_end,     'enable', 'off');
                
                % Rewind if necessary
                if tstep==2
                    tstep   = no_tstep;
                end
                % Plot
                for ii = tstep:-1:2
                    % Check for a status of an interrupt flag
                    flag = getappdata(sill1d_gui_handle, 'flag');
                    if flag == 1
                        setappdata(sill1d_gui_handle, 'flag', 0);
                        break;
                    end
                    
                    tic;
                    plot_results(ii);
                    drawnow;
                    % In case of fast update, pause for a second
                    timepstep = toc;
                    %if timepstep < .5
                    %    pause(0.5-timepstep);
                    %end
                end
                
                set(obj.sill1d_play_go_start,   'enable', 'on');
                set(obj.sill1d_play_backward,   'enable', 'on');
                set(obj.sill1d_play_stop,       'enable', 'on');
                set(obj.sill1d_play_forward,    'enable', 'on');
                set(obj.sill1d_play_go_end,     'enable', 'on');
                
            case 'sill1d_play_stop'
                % Set an interrupt flag and write into storage
                flag = 1;
                setappdata(sill1d_gui_handle, 'flag', flag);
                set(obj.sill1d_play_go_start,   'enable', 'on');
                set(obj.sill1d_play_backward,   'enable', 'on');
                set(obj.sill1d_play_stop,       'enable', 'on');
                set(obj.sill1d_play_forward,    'enable', 'on');
                set(obj.sill1d_play_go_end,     'enable', 'on');
                
            case 'sill1d_play_forward'
                % Set an interrupt flag and write into storage
                flag = 0;
                setappdata(sill1d_gui_handle, 'flag', flag);
                
                %  Enable appropreate play buttons
                set(obj.sill1d_play_go_start,   'enable', 'off');
                set(obj.sill1d_play_backward,   'enable', 'off');
                set(obj.sill1d_play_stop,       'enable', 'on' );
                set(obj.sill1d_play_forward,    'enable', 'off');
                set(obj.sill1d_play_go_end,     'enable', 'off');
                
                % Rewind if necessary
                if tstep==no_tstep
                    tstep   = 2;
                end
                % Plot
                for ii = tstep:no_tstep
                    % Check for a status of an interrupt flag
                    flag = getappdata(sill1d_gui_handle, 'flag');
                    if flag == 1
                        setappdata(sill1d_gui_handle, 'flag', 0);
                        break;
                    end
                    
                    tic;
                    plot_results(ii);
                    drawnow;
                    % In case of fast update, pause for a second
                    timepstep = toc;
                    %if timepstep < .5
                    %    pause(0.5-timepstep);
                    %end
                end
                
                set(obj.sill1d_play_go_start,   'enable', 'on');
                set(obj.sill1d_play_backward,   'enable', 'on');
                set(obj.sill1d_play_stop,       'enable', 'on');
                set(obj.sill1d_play_forward,    'enable', 'on');
                set(obj.sill1d_play_go_end,     'enable', 'on');
                
            case 'sill1d_play_go_end'
                plot_results(min(tstep+1, no_tstep));
                
            case 'sill1d_time'
                % Read current value
                time_ma     = str2double(get(obj.sill1d_time, 'String'));
                
                % Bail out if this does not evaluate to number
                if isnan(time_ma)
                    return;
                end
                
                % Find nearest timestep
                Time_comp   = (max(result.Time(:,1))-result.Time(:,1))/1e6;
                [~, ind]    = min(abs(Time_comp-time_ma));
                
                % Update plot
                if ~isempty(ind)
                    if ind==1
                        ind = 2;
                    end
                    plot_results(ind);
                end
        end
    end

%% fun plot_results
    function plot_results(tstep)
        % Store current tstep
        setappdata(sill1d_gui_handle, 'tstep', tstep);
        
        % Time in Ma
        time_ma     = round((max(result.Time(:,1))-result.Time(tstep,1))/1e6, 6);
        
        % Update GUI
        set(obj.sill1d_time,      'String', sprintf('%0.6f', time_ma)); % One year precision
        set(obj.sill1d_time_info, 'String', [' [Ma] (', num2str(tstep), '/', num2str(no_tstep), ')']);
        
        % Read Nodes that are active during this time and display those
        % only; required so that eroded layers do not show up in plots
        Active_nodes = find(result.Active(:,tstep)==1);
        
        % Finds steps where sills need to be plotted and when need to be removed
        ind_sill    = Active_nodes(result.Ind(Active_nodes)<0);
        no_sill     = unique(result.Ind(ind_sill));
        plots       = result;
        Ix          = [];
        for j = 1:length(no_sill)
            ind = find(result.Ind(Active_nodes)==no_sill(j));
            if length(ind)<2
                plots.Toc(Active_nodes(ind), tstep)         = NaN;
                plots.CO2_org(Active_nodes(ind), tstep)     = NaN;
                plots.CO2_release(Active_nodes(ind), tstep) = NaN;
            else
                Ix      = [Ix; no_sill(j)];
            end
        end
        
        %% - Temperature
        p_p1    = getappdata(sill1d_gui_handle, 'p_p1');
        if isempty(p_p1) || ~ishandle(p_p1)
            p_p1 = plot(result.Temp(Active_nodes, tstep), -result.Gcoord(Active_nodes, tstep)./1e3, '.-', 'Color', [1 0 0], 'Parent', ax.res(1)); % red
            add_menu(p_p1);
            ylabel(ax.res(1), 'TVDSS [km]');
            xlabel(ax.res(1), sprintf('Temperature [%cC]', char(176)));
            grid(  ax.res(1), 'on');
            set(   ax.res(1), 'xaxisLocation', 'top')
            setappdata(sill1d_gui_handle, 'p_p1', p_p1);
        else
            set(p_p1, 'xdata', result.Temp(Active_nodes, tstep), 'ydata', -result.Gcoord(Active_nodes, tstep)./1e3);
        end
        
        % Sets limits to remove plotting effects of minor inaacuracies when constant T is used
        xlim(  ax.res(1), [min(p_p1.XData)-1 max(p_p1.XData)+1]);
        if abs(max(p_p1.XData)-min(p_p1.XData))<2
            xlim(ax.res(1), [min(p_p1.XData)-2 max(p_p1.XData)+2]);
        else
            xlim(  ax.res(1), [min(p_p1.XData)-1 max(p_p1.XData)+1]);
            if min(p_p1.XData)>=-1e-2
                xlims = xlim(ax.res(1));
                xlims(1) = 0;
                xlim(ax.res(1), xlims);
            end
        end
        ylim(  ax.res(1), [min(-result.Gcoord(Active_nodes, tstep)./1e3) max(-result.Gcoord(Active_nodes, tstep)./1e3)]);       
        
        % Plots Sills if present
        hold(ax.res(1), 'on')
        delete(findobj(ax.res(1), 'Type', 'Patch')); %Sills are plotted everytime and have to be manually removed
        for k = 1:length(Ix)
            X = [ax.res(1).XLim(1) ax.res(1).XLim(2) ax.res(1).XLim(2) ax.res(1).XLim(1)];
            Y = [max(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 max(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 min(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 min(result.Gcoord(result.Ind==Ix(k), tstep))/1e3];
            patch(X, -Y, [224/255 224/255 224/255], 'FaceAlpha', .5, 'EdgeColor', 'none', 'Parent', ax.res(1));
        end
        hold(ax.res(1), 'off')
        
        %% - Vitrinite
        p_p2    = getappdata(sill1d_gui_handle, 'p_p2');
        if isempty(p_p2) || ~ishandle(p_p2)
            p_p2 = plot(result.Ro(Active_nodes, tstep), -result.Gcoord(Active_nodes, tstep)./1e3, '.-', 'Color', [127/255 0 255/255], 'Parent', ax.res(2)); % purple
            add_menu(p_p2);
            ylabel(ax.res(2), 'TVDSS [km]')
            xlabel(ax.res(2), 'Vitrinite Reflectance [%Ro]');
            set(   ax.res(2), 'xaxisLocation', 'top')
            grid(  ax.res(2), 'on')
            setappdata(sill1d_gui_handle, 'p_p2', p_p2);
        else
            set(p_p2, 'xdata', result.Ro(Active_nodes, tstep), 'ydata', -result.Gcoord(Active_nodes, tstep)./1e3);
        end
        ylim(  ax.res(2), [min(-result.Gcoord(Active_nodes, tstep)./1e3) max(-result.Gcoord(Active_nodes, tstep)./1e3)]);
        xlim(  ax.res(2), 'auto');
               
        % Plots VR data if available for last step
        delete(findobj(ax.res(2), 'Color', 'k')); % Removes plotted VR data
        delete(findobj(ax.res(2), 'Color', 'b')); % Removes plotted VR data
        if tstep == size(result.Temp, 2) && ~isempty(welldata.VR)
            hold(ax.res(2), 'on');
            if size(welldata.VR,2)<3
               plot(welldata.VR(:,2), -welldata.VR(:,1)./1e3, 'o', 'Color', [0 0 0], 'Parent', ax.res(2));
            else
               errorbarxy(ax.res(2), welldata.VR(:,2), -welldata.VR(:,1)./1e3, welldata.VR(:,3),0.*welldata.VR(:,3),{'ko', 'b', 'b'});
            end
            hold(ax.res(2), 'off');
            xlim(ax.res(2), [min([welldata.VR(:,2); result.Ro(Active_nodes, tstep)]) max([welldata.VR(:,2); result.Ro(Active_nodes, tstep)])]);
        end
         
        % Plots Sills if present
        hold(ax.res(2), 'on')
        delete(findobj(ax.res(2), 'Type', 'Patch')); %Sills are plotted everytime and have to be manually removed
        for k = 1:length(Ix)
            X = [ax.res(2).XLim(1) ax.res(2).XLim(2) ax.res(2).XLim(2) ax.res(2).XLim(1)];
            Y = [max(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 max(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 min(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 min(result.Gcoord(result.Ind==Ix(k), tstep))/1e3];
            patch(X, -Y, [224/255 224/255 224/255], 'FaceAlpha',.5, 'EdgeColor', 'none', 'Parent', ax.res(2));
            xlim(ax.res(2), X(1:2));
        end
        hold(ax.res(2), 'off')
        
        %% - TOC
        p_p3    = getappdata(sill1d_gui_handle, 'p_p3');
        if isempty(p_p3) || ~ishandle(p_p3)
            p_p3 = plot(plots.Toc(Active_nodes,tstep), -result.Gcoord(Active_nodes, tstep)./1e3, '.-', 'Color', [0 153/255 0], 'Parent', ax.res(3)); % Green
            add_menu(p_p3);
            ylabel(ax.res(3), 'TVDSS [km]');
            xlabel(ax.res(3), 'TOC Content [wt%]');
            set(   ax.res(3), 'xaxisLocation', 'top')
            grid(  ax.res(3), 'on');
            setappdata(sill1d_gui_handle, 'p_p3', p_p3);
        else
            set(p_p3, 'xdata', plots.Toc(Active_nodes,tstep), 'ydata', -result.Gcoord(Active_nodes, tstep)./1e3);
        end
        ylim(  ax.res(3), [min(-result.Gcoord(Active_nodes, tstep)./1e3) max(-result.Gcoord(Active_nodes, tstep)./1e3)]);
        xlim(  ax.res(3), 'auto');

        % Plots TOC data if available for last step
        delete(findobj(ax.res(3), 'Color', 'k')); % Removes plotted TOC data
        delete(findobj(ax.res(3), 'Color', 'b')); % Removes plotted TOC data
        if tstep == size(result.Temp, 2) && ~isempty(welldata.Toc)
            hold(ax.res(3), 'on');
            if size(welldata.Toc,2)<3
                h_p = plot(welldata.Toc(:,2), -welldata.Toc(:,1)./1e3, 'o', 'Color', [0 0 0], 'Parent', ax.res(3));
            else
                h_p = errorbarxy(ax.res(2), welldata.Toc(:,2), -welldata.Toc(:,1)./1e3, welldata.Toc(:,3),0.*welldata.Toc(:,3),{'ko', 'b', 'b'});
            end
            hold(ax.res(3), 'off');
            xlim(ax.res(3), [min([welldata.Toc(:,2); plots.Toc(Active_nodes,tstep)]) max([welldata.Toc(:,2); plots.Toc(Active_nodes,tstep)])]);
        end
        
        % Plots Sills if present
        hold(ax.res(3), 'on')
        delete(findobj(ax.res(3), 'Type', 'Patch')); %Sills are plotted everytime and have to be manually removed
        for k = 1:length(Ix)
            X = [ax.res(3).XLim(1) ax.res(3).XLim(2) ax.res(3).XLim(2) ax.res(3).XLim(1)];
            Y = [max(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 max(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 min(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 min(result.Gcoord(result.Ind==Ix(k), tstep))/1e3];
            patch(X, -Y, [224/255 224/255 224/255], 'FaceAlpha', .5, 'EdgeColor', 'none', 'Parent', ax.res(3));
        end
        hold(ax.res(3), 'off')
               
        %% T-Max
        p_p4    = getappdata(sill1d_gui_handle, 'p_p4');
        T_max               = result.Tmax(:,tstep);
        T_max(result.Ind<0) = NaN;
        if isempty(p_p4) || ~ishandle(p_p4)
            p_p4 = plot(T_max(Active_nodes), -result.Gcoord(Active_nodes, tstep)./1e3, '.-', 'Color', [255/255 128/255 0], 'MarkerFaceColor', 'k', 'Parent', ax.res(4)); % Orange
            add_menu(p_p4);
            ylabel(ax.res(4), 'TVDSS [km]');
            xlabel(ax.res(4), sprintf('T_m_a_x [%cC]', char(176)));
            set(   ax.res(4), 'xaxisLocation', 'top')
            grid(  ax.res(4), 'on');
            setappdata(sill1d_gui_handle, 'p_p4', p_p4);
        else
            set(p_p4, 'xdata', T_max(Active_nodes), 'ydata', -result.Gcoord(Active_nodes, tstep)./1e3);
        end
        
        % Sets limits to remove plotting effects of minor inaacuracies when constant T is used
        if abs(max(p_p4.XData)-min(p_p4.XData))<2
            xlim(ax.res(4), [min(p_p4.XData)-2 max(p_p4.XData)+2]);
        else
            xlim(ax.res(4), 'auto');
        end
        ylim(  ax.res(4), [min(-result.Gcoord(Active_nodes, tstep)./1e3) max(-result.Gcoord(Active_nodes, tstep)./1e3)]);
        
        % Plots Sills if present
        hold(ax.res(4), 'on')
        delete(findobj(ax.res(4), 'Type', 'Patch')); %Sills are plotted everytime and have to be manually removed
        for k = 1:length(Ix)
            X = [ax.res(4).XLim(1) ax.res(4).XLim(2) ax.res(4).XLim(2) ax.res(4).XLim(1)];
            Y = [max(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 max(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 min(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 min(result.Gcoord(result.Ind==Ix(k), tstep))/1e3];
            patch(X, -Y, [224/255 224/255 224/255], 'FaceAlpha', .5, 'EdgeColor', 'none', 'Parent', ax.res(4));
        end
        hold(ax.res(4), 'off')
        
        %% Pressure
        p_p5    = getappdata(sill1d_gui_handle, 'p_p5');
        if isempty(p_p5) || ~ishandle(p_p5)
            p_p5 = plot(result.Pres(Active_nodes,tstep)./1e6, -result.Gcoord(Active_nodes,tstep)./1e3, '.-', 'Color', [0 24/255 204/255], 'Parent', ax.res(5)); % Blue
            add_menu(p_p5);
            ylabel(ax.res(5), 'TVDSS [km]')
            xlabel(ax.res(5), 'Pressure [MPa]');
            set(   ax.res(5), 'xaxisLocation', 'top')
            grid(  ax.res(5), 'on')
            setappdata(sill1d_gui_handle, 'p_p5', p_p5);
        else
            set(p_p5, 'xdata', result.Pres(Active_nodes,tstep)./1e6, 'ydata', -result.Gcoord(Active_nodes, tstep)./1e3);
        end
        ylim(  ax.res(5), [min(-result.Gcoord(Active_nodes, tstep)./1e3) max(-result.Gcoord(Active_nodes, tstep)./1e3)])

        % Plots Sills if present
        hold(ax.res(5), 'on')
        delete(findobj(ax.res(5), 'Type', 'Patch')); %Sills are plotted everytime and have to be manually removed
        for k = 1:length(Ix)
            X = [ax.res(5).XLim(1) ax.res(5).XLim(2) ax.res(5).XLim(2) ax.res(5).XLim(1)];
            Y = [max(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 max(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 min(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 min(result.Gcoord(result.Ind==Ix(k), tstep))/1e3];
            patch(X, -Y, [224/255 224/255 224/255], 'FaceAlpha', .5, 'EdgeColor', 'none', 'Parent', ax.res(5));
        end
        hold(ax.res(5), 'off')
        
        %% - CO2 Inorganic
        p_p6    = getappdata(sill1d_gui_handle, 'p_p6');
        if isempty(p_p6) || ~ishandle(p_p6)
            p_p6 = plot(plots.CO2_release(Active_nodes,tstep), -result.Gcoord(Active_nodes,tstep)./1e3, '.-', 'Color', [253/255 153/255 204/255], 'Parent', ax.res(6)); % Pink
            add_menu(p_p6);
            ylabel(ax.res(6), 'TVDSS [km]');
            xlabel(ax.res(6), 'CO_2 inorganic [kg/m^3]');
            set(   ax.res(6), 'xaxisLocation', 'top')
            grid(  ax.res(6), 'on');
            setappdata(sill1d_gui_handle, 'p_p6', p_p6);
        else
            set(p_p6, 'xdata', plots.CO2_release(Active_nodes,tstep), 'ydata', -result.Gcoord(Active_nodes, tstep)./1e3);
        end
        ylim(  ax.res(6), [min(-result.Gcoord(Active_nodes, tstep)./1e3) max(-result.Gcoord(Active_nodes, tstep)./1e3)]);

        % Plots Sills if present
        hold(ax.res(6), 'on')
        delete(findobj(ax.res(6), 'Type', 'Patch')); %Sills are plotted everytime and have to be manually removed
        for k = 1:length(Ix)
            X = [ax.res(6).XLim(1) ax.res(6).XLim(2) ax.res(6).XLim(2) ax.res(6).XLim(1)];
            Y = [max(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 max(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 min(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 min(result.Gcoord(result.Ind==Ix(k), tstep))/1e3];
            patch(X, -Y, [224/255 224/255 224/255], 'FaceAlpha', .5, 'EdgeColor', 'none', 'Parent', ax.res(6));
        end
        hold(ax.res(6), 'off')
        
        %% - CO2 Organic
        p_p7    = getappdata(sill1d_gui_handle, 'p_p7');
        if isempty(p_p7) || ~ishandle(p_p7)
            p_p7 = plot(plots.CO2_org(Active_nodes,tstep), -result.Gcoord(Active_nodes,tstep)./1e3, '.-', 'Color', [0 0 0], 'Parent', ax.res(7)); % Pink
            add_menu(p_p7);
            ylabel(ax.res(7), 'TVDSS [km]');
            xlabel(ax.res(7), 'CO_2 Organic [kg/m^3]');
            set(   ax.res(7), 'xaxisLocation', 'top')
            grid(  ax.res(7), 'on');
            setappdata(sill1d_gui_handle, 'p_p7', p_p7);
        else
            set(p_p7, 'xdata', plots.CO2_org(Active_nodes,tstep), 'ydata', -result.Gcoord(Active_nodes, tstep)./1e3);
        end
        ylim(  ax.res(7), [min(-result.Gcoord(Active_nodes, tstep)./1e3) max(-result.Gcoord(Active_nodes, tstep)./1e3)]);

        % Plots Sills if present
        hold(ax.res(7), 'on')
        delete(findobj(ax.res(7), 'Type', 'Patch')); %Sills are plotted everytime and have to be manually removed
        for k = 1:length(Ix)
            X = [ax.res(7).XLim(1) ax.res(7).XLim(2) ax.res(7).XLim(2) ax.res(7).XLim(1)];
            Y = [max(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 max(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 min(result.Gcoord(result.Ind==Ix(k), tstep))/1e3 min(result.Gcoord(result.Ind==Ix(k), tstep))/1e3];
            patch(X, -Y, [224/255 224/255 224/255], 'FaceAlpha', .5, 'EdgeColor', 'none', 'Parent', ax.res(7));
        end
        hold(ax.res(7), 'off')
    end

%% fun save_hires_fig
    function save_hires_fig()
        
        %  Open dialogbox
        [filename, pathname] = uiputfile(...
            {'*.png'},...
            'Save as');
        
        % If any folder is selected
        if length(filename)>1
            % Figure out which tab is active
            switch tab_panel.Selection
                case 1
                    h_active    = h_tp_input;
                case 2
                    h_active    = h_t2_hb2;
                case 3
                    h_active    = h_tp_release;
            end
            
            Position    = getpixelposition(h_active);
            
            % Colormap
            Map         = colormap;
            
            % Make hidden figure of the same size
            sill1d_hidden_plot = figure('units','pixels', 'pos', Position,...
                'NumberTitle','off', 'Name','Sill1D Hidden Plot','tag','sill1d_hidden_plot',...
                'MenuBar','none',...
                'PaperPositionMode','auto', ...
                'visible', 'off', ...
                'Color', [1 1 1]);
                               
            % Copy the container to the figure
            new_handle  = copyobj(h_active, sill1d_hidden_plot);
            
            % Activate colormap
            colormap(sill1d_hidden_plot, Map);
            
            % White background on all uicontainers
            h_uc    = findall(sill1d_hidden_plot, 'type', 'UIContainer');
            for i=1:length(h_uc)
                h_uc(i).BackgroundColor     = 'w';
            end
            
            % Print it
            set(sill1d_gui_handle, 'pointer', 'watch');
            try
                print(sill1d_hidden_plot, '-dpng', '-r150', [pathname, filename]);
            catch
                % Oh well, at least we can stop the spinner
                disp(' Figure save failed!');
            end
            set(sill1d_gui_handle, 'pointer', 'arrow');
            
            % Delete figure
            delete(sill1d_hidden_plot);            
        end
        
        disp(' Figure saved');
    end
end

%% fun data2clipboard
function data2clipboard(~, ~, obj)
DATA    = [obj.XData', obj.YData'];
num2clip(DATA);
end

%% fun add_menu
function add_menu(object_handle)
uicm    = uicontextmenu;
uimenu(uicm, 'Label', 'Copy data to clipboard', 'Callback', {@data2clipboard, object_handle});
object_handle.UIContextMenu = uicm;
end