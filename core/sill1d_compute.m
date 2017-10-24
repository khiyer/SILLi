function [Results, Time] = sill1d_compute(ROCK, SILL, FLUID, WellData, Res, codeloc)
% sill1d_comput
%
% 1D Thermal cooling of intrusive
%
% Developed by Karthik Iyer, Henrik Svensen and Daniel W. Schmid
%
%--------------------------------------------------------------------------

%% Initialization
disp(' - Setting up numerical model');
yr              = 365*24*60*60;

%% - Generate Full Mesh for all Lithologies
Gcoord          = [];

% Mesh generated based on user defined resolution
basement        = max([ROCK.Tops+10; SILL.Tops+SILL.Thick+(SILL.Thick.*5)]);
All_tops        = unique([ROCK.Tops; basement; SILL.Tops; SILL.Tops+SILL.Thick]);
for i = 1:length(All_tops)-1
    if any(All_tops(i)==SILL.Tops);
        d_z         = min([Res.dz_sill (All_tops(i+1)-All_tops(i))/Res.pts_sill]);
        Gcoord      = unique([Gcoord All_tops(i):d_z:All_tops(i+1) All_tops(i+1)]);
    else
        d_z         = min([Res.dz_sed (All_tops(i+1)-All_tops(i))/Res.pts_sed]);
        Gcoord      = unique([Gcoord All_tops(i):d_z:All_tops(i+1) All_tops(i+1)]);
    end
end

nnod            = length(Gcoord);
nel             = nnod-1;
ELEM2NODE       = [1:nnod-1; 2:nnod];

% Generate Indices for different sediments
Ind             = NaN*ones(length(Gcoord),1);
Ind_nel         = NaN*ones(nel,1);
Gcoord_c        = (Gcoord(1:end-1)+Gcoord(2:end))./2;
for i=1:length(ROCK.Tops)
    Ind(Gcoord>=ROCK.Tops(i)) = i;
    Ind_nel(Gcoord_c>=ROCK.Tops(i)) = i;
end

% Generate Indices for sills
sill_ind        = 1;
Events_sill     = sort(unique([SILL.E_time]), 'descend');
for i=1:length(Events_sill)
    x           = find(SILL.E_time==Events_sill(i));
    for j = 1:length(x)
        Ind(Gcoord>=SILL.Tops(x(j)) & Gcoord<(SILL.Tops(x(j))+SILL.Thick(x(j)))) = -sill_ind;
        Ind_nel(Gcoord_c>=SILL.Tops(x(j)) & Gcoord_c<(SILL.Tops(x(j))+SILL.Thick(x(j)))) = -sill_ind;
        sill_ind = sill_ind+1;
    end
end

%% - Generate Time Sequence for element depostion and sill emplacement

% Takes the top and bottom ages of a layer together with resolution and creates
% deposition time for each node
Events_dep      = sort(unique([ROCK.Ages; ROCK.Ages(end)+1e-6]), 'descend');
Times           = zeros(nel,1);
water_depth     = ROCK.Tops(1);
Gcoord_n        = water_depth*ones(nnod,1);
Dx              = zeros(nel,1);
Ero_times       = ROCK.Ero_t(~isnan(ROCK.Ero_t));
for i = 2:length(Events_dep)
    x                       = find(ROCK.Ages==Events_dep(i));
    Selected                = Gcoord(ELEM2NODE(:,Ind_nel==x));
    dx                      = Selected(2,:)-Selected(1,:);
    Dx(Ind_nel==x)          = dx';
    Gcoord_sel              = [0 cumsum(dx)];
    ero_check               = find(Ero_times<Events_dep(i-1) & Ero_times>Events_dep(i));
    if isempty(ero_check)
        t                       = interp1([max(Gcoord_sel) min(Gcoord_sel)], [Events_dep(i-1) Events_dep(i)], Gcoord_sel);
        Times(Ind_nel==x)       = t(1:end-1)';
    else
        t                       = interp1([max(Gcoord_sel) min(Gcoord_sel)], [min(Ero_times(ero_check)) Events_dep(i)], Gcoord_sel);
        Times(Ind_nel==x)       = t(1:end-1)';
    end
end

% Sets sill intrusion times for all nodes within sills
sill_ind        = 1;
for i = 1:length(Events_sill)
    x                       = find(SILL.E_time==Events_sill(i));
    for j = 1:length(x)
        Times(Ind_nel==-sill_ind)   = Events_sill(i);
        Selected                    = Gcoord(ELEM2NODE(:,Ind_nel==-sill_ind));
        dx                          = Selected(2,:)-Selected(1,:);
        Dx(Ind_nel==-sill_ind)      = dx';
        sill_ind                    = sill_ind+1;
    end
end

%% - Check if sediment deposition time exactly matches sill emplacement time
for i = 1:length(Events_sill)
    ind_dep     = find(Times==Events_sill(i) & Ind_nel>0);
    if ~isempty(ind_dep)
        Times(ind_dep) = (Events_sill(i)+Times(ind_dep+1))/2;
    end
end

%% - Check if sediment deposition time exactly matches erosion time
Times_ero       = unique(ROCK.Ero_t(ROCK.Ero_t>0));
for i = 1:length(Times_ero)
    ind_dep     = find(Times==Times_ero(i) & Ind_nel>0);
    if ~isempty(ind_dep)
        Times(ind_dep) = (Times_ero(i)+Times(ind_dep+1))/2;
    end
end

%% - Construct Time Sequence
Times_all       = sort(([unique(Times); unique(ROCK.Ero_t(ROCK.Ero_t>0))]), 'descend');
if min(Times_all)~=0
    if ismember(min(Times_all), SILL.E_time) || ismember(min(Times_all), Times_ero)
        Times_all       = [Times_all; 0];
    else
        Times_all       = [Times_all; min(Times_all)-1e-8; 0];
    end
else
    add_time        = min([min(Times(Times>0))/2 1e-8]);
    Times_all       = sort([Times_all; add_time], 'descend');
    Times(1)        = add_time;
end
time_all        = 0;

%% - Initialize Rock Properties
R                   = 8.31;
A                   = 1e13;
Ro                  = NaN*ones(nnod,1);
Rho_rock            = NaN*ones(nnod,1);
Cp_rock             = NaN*ones(nnod,1);
Rho_melt            = NaN*ones(nnod,1);
Cp_melt             = NaN*ones(nnod,1);
Phi_rock            = NaN*ones(nnod,1);
K_rock              = NaN*ones(nnod,1);
Toc_rock            = NaN*ones(nnod,1);
Toc_rock_o          = NaN*ones(nnod,1);
Lm_rock             = NaN*ones(nnod,1);
Ld_rock             = NaN*ones(nnod,1);
Lc_rock             = NaN*ones(nnod,1);
Sol_rock            = NaN*ones(nnod,1);
Liq_rock            = NaN*ones(nnod,1);
Latent_dehyd        = NaN*ones(nnod,1);
Latent_om           = NaN*ones(nnod,1);
W                   = NaN*ones(nnod,20);
for i = 1:ROCK.num
    Rho_rock(Ind==i)    = ROCK.Rho(i);
    Cp_rock(Ind==i)     = ROCK.Cp(i);
    Phi_rock(Ind==i)    = ROCK.Phi(i);
    K_rock(Ind==i)      = ROCK.K(i);
    Toc_rock(Ind==i)    = ROCK.Toc(i);
    Toc_rock_o(Ind==i)  = ROCK.Toc(i);
    Lm_rock(Ind==i)     = ROCK.Lm(i);
    Ld_rock(Ind==i)     = ROCK.Ld(i);
    Lc_rock(Ind==i)     = 0;
    Sol_rock(Ind==i)    = 0;
    Liq_rock(Ind==i)    = 0;
    Latent_dehyd(Ind==i)= 0;
    Latent_om(Ind==i)   = 0;
    W(Ind==i,:)         = 1;
end
for i = 1:SILL.num
    Rho_rock(Ind==-i)    = SILL.Rhos(i);
    Cp_rock(Ind==-i)     = SILL.Cps(i);
    Rho_melt(Ind==-i)    = SILL.Rhom(i);
    Cp_melt(Ind==-i)     = SILL.Cpm(i);
    Phi_rock(Ind==-i)    = 0;
    K_rock(Ind==-i)      = SILL.K(i);
    Toc_rock(Ind==-i)    = 0;
    Toc_rock_o(Ind==-i)  = 0;
    Lm_rock(Ind==-i)     = 0;
    Ld_rock(Ind==-i)     = 0;
    Lc_rock(Ind==-i)     = SILL.Lc(i);
    Sol_rock(Ind==-i)    = SILL.Sol(i);
    Liq_rock(Ind==-i)    = SILL.Liq(i);
    Latent_dehyd(Ind==-i)= 0;
    Latent_om(Ind==-i)   = 0;
    W(Ind==-i,:)         = 0;
end

CO2_org         = zeros(nnod,1);
K_solid         = (K_rock.^(1-Phi_rock)).*(FLUID.k.^(Phi_rock));

%% - Load Carbonate Data
CO2_release     = zeros(nnod,1);
Ind_carb        = zeros(nnod,1);
if max(ROCK.Carb)>0
    load(fullfile(codeloc, 'dat', 'Marl.mat'))
    load(fullfile(codeloc, 'dat', 'Dolostone.mat'))
    load(fullfile(codeloc, 'dat', 'DolostoneEvaporite.mat'))
    
    C1_lith     = find(ROCK.Carb==1);
    for a = 1:length(C1_lith)
        Ind_carb(Ind==C1_lith(a))   = 1;
    end
    C2_lith     = find(ROCK.Carb==2);
    for a = 1:length(C2_lith)
        Ind_carb(Ind==C2_lith(a))   = 2;
    end
    C3_lith     = find(ROCK.Carb==3);
    for a = 1:length(C3_lith)
        Ind_carb(Ind==C3_lith(a))   = 3;
    end
    
    C1                      = scatteredInterpolant(Marl.T(:), Marl.P(:).*1e5, Marl.CO2(:));
    C2                      = scatteredInterpolant(Dolo.T(:), Dolo.P(:).*1e5, Dolo.CO2(:));
    C3                      = scatteredInterpolant(Dol_ev.T(:), Dol_ev.P(:).*1e5, Dol_ev.CO2(:));
end

%% - Initialize Temperature
Temp            = WellData.T(1,2)*ones(nnod,1);
T_check         = NaN*ones(size(Temp));
T_max           = Temp;

x               = WellData.T(:,1)';
y               = WellData.T(:,2)';
polya           = polyfit(x,y,1);

%% - Initialize Pressure
Pres            = zeros(nnod,1);

%% - Initialize Output Variables
Results.nel             = nel;          % No. of elements
Results.nnod            = nnod;         % No. of nodes
Results.Gcoord_c        = Gcoord_c;     % Z-values of element centres
Results.Ind             = Ind;          % Lithology indices (node-based)
Results.Ind_nel         = Ind_nel;      % Lithology indices (element-based)
Results.Ind_carb        = Ind_carb;     % Carbonate indices (node-based)

Results.Gcoord(:,1)     = Gcoord;       % Z-values of entire column including eroded layers
Results.Temp(:,1)       = Temp;         % Temperature (?C)
Results.Pres(:,1)       = Pres;         % Pressure (Pa)
Results.Toc(:,1)        = Toc_rock.*1e2;% TOC contents (wt%)
Results.CO2_org(:,1)    = CO2_org;      % Organic CO2 generated (kg/m3)
Results.Ro(:,1)         = Ro;           % Maturity (%Ro)
Results.Tmax(:,1)       = T_max;        % Maximum Temperature (?C)
Results.Active(:,1)     = zeros(nnod,1);% Active nodes for time step
Results.CO2_release(:,1)= CO2_release;  % CO2 release (kg/m3)
Results.Time(1,1)       = 0;           % Time (yr)
c1                      = 2;

Active          = false(nel,1);

%% Time Loop
textprogressbar(' - Time Loop: ');
for i = 1:length(Times_all)-1
    textprogressbar(i/length(Times_all)*100);
    
    s_emp   = find(unique(SILL.E_time==Times_all(i)));
    ero     = find(unique(ROCK.Ero_t(~isnan(ROCK.Ero_t)))==Times_all(i));
    
    if isempty(s_emp) && isempty(ero)
        %%- Deposition
        
        % Set active layer that is deposited
        active_now                  = find(Times==Times_all(i));
        Gcoord_n(active_now+1:end)  = Gcoord_n(active_now+1:end) + Dx(active_now);
        Active(Times==Times_all(i)) = 1;
        
        % Create new mesh only for thermal diffusion
        n_nel = size(ELEM2NODE(:,Active'),2);
        E2N   = ELEM2NODE(:,Active');
        L2G = [1:n_nel; 2:n_nel+1];
        Active_nodes = unique(E2N(:));
        
        Pres(Active_nodes(1))       = min(Gcoord).*FLUID.rho*9.81 + 1e5;
        Pres(Active_nodes(2:end))   = cumsum(Rho_rock(Active_nodes(1:end-1)).*diff(Gcoord_n(Active_nodes)).*9.81) + Pres(Active_nodes(1));
        
        T_top                   = WellData.T(1,2);
        T_bot                   = (max(Gcoord_n(E2N(:)))).*polya(1)+polya(2);
        
        % Set time allowed for VR calculation
        if i < length(Times_all)
            dt_vr         = (Times_all(i)-Times_all(i+1))*1e6*yr;
        else
            dt_vr         = 1e-2;
        end
        
        time_all = time_all + dt_vr/yr;
        
        % Recompute Densities and heat capacities for coefficient in thermal
        % solver
        Rho_solid               = (1-Phi_rock).*Rho_rock + Phi_rock.*FLUID.rho;
        Cp_solid                = (1-Phi_rock).*Cp_rock + Phi_rock.*FLUID.cp;
        
        % Change heat capacity if sill T is between liquidus and solidus
        solid_ind               = find(Ind<0 & Temp<=Sol_rock);
        Cp_solid(solid_ind)     = Cp_rock(solid_ind);
        melt_solid              = find(Ind<0 & Temp>Sol_rock & Temp<Liq_rock);
        Rho_solid(melt_solid)   = Rho_melt(melt_solid);
        Cp_solid(melt_solid)    = Cp_melt(melt_solid).*(1 + (Lc_rock(melt_solid)./((Liq_rock(melt_solid)-Sol_rock(melt_solid)).*Cp_melt(melt_solid))));
        melt_ind                = find(Ind<0 & Temp>=Liq_rock);
        Rho_solid(melt_ind)     = Rho_melt(melt_ind);
        Cp_solid(melt_ind)      = Cp_melt(melt_ind);
        
        Coeff_temp              = Rho_solid.*Cp_solid;
        
        
        Temp_old = Temp;
        T_max(Temp>T_max) = Temp(Temp>T_max);
        
        % Solve for Diffusion
        T = thermal1d_fem(Temp, Gcoord_n, Coeff_temp, K_solid, T_top, T_bot, dt_vr, n_nel, n_nel+1, Latent_dehyd, Latent_om,...
            L2G, E2N);
        
        % Reorder Temperature in Global Mesh
        T_check(E2N(:)) = T(L2G(:));
        Temp(E2N(:)) = T(L2G(:));
        
        T_max(Temp>T_max) = Temp(Temp>T_max);
        
        %Calculate CO2 release from carbonates
        if max(ROCK.Carb)>0
            C1_ind      = find(Ind_carb(Active_nodes)==1);
            C2_ind      = find(Ind_carb(Active_nodes)==2);
            C3_ind      = find(Ind_carb(Active_nodes)==3);
            
            if ~isempty(C1_ind)
                C1_now  = (1-Phi_rock(Active_nodes(C1_ind))).*Rho_rock(Active_nodes(C1_ind)).*(C1(Temp(Active_nodes(C1_ind)), Pres(Active_nodes(C1_ind)))./100);
                C1_do  = find(CO2_release(Active_nodes(C1_ind))<C1_now);
                CO2_release(Active_nodes(C1_ind(C1_do))) = C1_now(C1_do);
            end
            if ~isempty(C2_ind)
                C2_now = (1-Phi_rock(Active_nodes(C2_ind))).*Rho_rock(Active_nodes(C2_ind)).*(C2(Temp(Active_nodes(C2_ind)), Pres(Active_nodes(C2_ind)))./100);
                C2_do  = find(CO2_release(Active_nodes(C2_ind))<C2_now);
                CO2_release(Active_nodes(C2_ind(C2_do))) = C2_now(C2_do);
            end
            if ~isempty(C3_ind)
                C3_now = (1-Phi_rock(Active_nodes(C3_ind))).*Rho_rock(Active_nodes(C3_ind)).*(C3(Temp(Active_nodes(C3_ind)), Pres(Active_nodes(C3_ind)))./100);
                C3_do  = find(CO2_release(Active_nodes(C3_ind))<C3_now);
                CO2_release(Active_nodes(C3_ind(C3_do))) = C3_now(C3_do);
            end
        end
        
        % Calculate maturation
        [Latent_dehyd(Active_nodes), Latent_om(Active_nodes), Toc_rock_o(Active_nodes), W(Active_nodes,:), Ro(Active_nodes)] = coeff_diff(Temp(Active_nodes),...
            dt_vr, Rho_rock(Active_nodes), A, R, length(Active_nodes), Ind(Active_nodes), Toc_rock(Active_nodes),Lm_rock(Active_nodes), Temp_old(Active_nodes),...
            Toc_rock_o(Active_nodes), Phi_rock(Active_nodes), W(Active_nodes,:), Ld_rock(Active_nodes));
        
        CO2_org(Active_nodes)     = (Toc_rock(Active_nodes)- Toc_rock_o(Active_nodes)).*(1-Phi_rock(Active_nodes)).*Rho_rock(Active_nodes).*3.66;
        
        Ro_plot = Ro;
        Ro_plot(Ind<0) = NaN;
        
        A_store                   = zeros(nnod,1);
        A_store(Active_nodes)     = 1;
        
        % Store Dynamic Variables
        Results.Temp(:,c1)            = Temp;                % Temperature
        Results.Pres(:,c1)            = Pres;                % Pressure
        Results.Toc(:,c1)             = Toc_rock_o.*1e2;     % TOC contents
        Results.Ro(:,c1)              = Ro_plot;             % Maturity
        Results.Tmax(:,c1)            = T_max;               % Maximum Temperature
        Results.Active(:,c1)          = A_store;             % Active nodes for time step
        Results.Gcoord(:,c1)          = Gcoord_n;            % Z-values for lithologies
        Results.CO2_release(:,c1)     = CO2_release;         % CO2 release
        Results.Time(c1,1)            = time_all;            % Time
        Results.CO2_org(:,c1)         = CO2_org;             % TOC reacted to CO2
        
        c1 = c1 + 1;
        
    elseif ~isempty(s_emp)
        
        %% - Sill emplacement
        a       = unique(Ind_nel(Times==Times_all(i)));
        b       = sort(a,'descend');
        bb      = SILL.E_temp(SILL.E_time==Times_all(i));
        sill_t  = max(SILL.Thick(SILL.E_time==Times_all(i)))/2;
        sill_d  = SILL.Tops(SILL.E_time==Times_all(i));
        rock_k  = ROCK.K(ROCK.Tops==max(ROCK.Tops(ROCK.Tops<sill_d(1))));
        rock_r  = ROCK.Rho(ROCK.Tops==max(ROCK.Tops(ROCK.Tops<sill_d(1))));
        rock_cp = ROCK.Cp(ROCK.Tops==max(ROCK.Tops(ROCK.Tops<sill_d(1))));
        rock_phi= ROCK.Phi(ROCK.Tops==max(ROCK.Tops(ROCK.Tops<sill_d(1))));
        rhocp   = (1-rock_phi)*rock_r*rock_cp + rock_phi*FLUID.rho*FLUID.cp;
        
        for l = 1:length(b)
            Temp(Ind==b(l)) = bb(l);
        end
        
        % Set active layer that is deposited
        for x = 1:length(b)
            active_now   = find(Ind_nel==b(x));
            for j= active_now(1)+1:active_now(end)
                Gcoord_n(j)     = Gcoord_n(j) + sum(Dx(active_now(1):j-1));
            end
            Gcoord_n(active_now(end)+1:end) = Gcoord_n(active_now(end)+1:end) + sum(Dx(active_now(1:end)));
        end
        
        Active(Times==Times_all(i)) = 1;
        
        % Create new mesh only for thermal diffusion
        n_nel = size(ELEM2NODE(:,Active'),2);
        E2N   = ELEM2NODE(:,Active');
        L2G = [1:n_nel; 2:n_nel+1];
        Active_nodes = unique(E2N(:));
        
        Pres(Active_nodes(1))       = min(Gcoord).*1.05e3*9.81;
        Pres(Active_nodes(2:end))   = cumsum(Rho_rock(Active_nodes(1:end-1)).*diff(Gcoord_n(Active_nodes)).*9.81) + Pres(Active_nodes(1));
        
        T_top                   = WellData.T(1,2);
        T_bot                   = (max(Gcoord_n(E2N(:)))).*polya(1)+polya(2);

        if i < length(Times_all)
            dt_vr         = (Times_all(i)-Times_all(i+1))*1e6;
        else
            dt_vr = 1e5;
        end
        time = 0;
        sill_time       = ceil(10*sill_t^2*rhocp/rock_k/yr); % t = 10tau*d^2/D;
        time_end_s      = min([sill_time dt_vr]);
        dt_range        = [time_end_s/1e5 time_end_s/10];
        time_range      = linspace(0,time_end_s,1000);
        a_exp           = dt_range(1);
        b_exp           = log(dt_range(2)/dt_range(1))/time_end_s;
        dt_s            = a_exp.*exp(b_exp.*time_range);
        dt_sill         = dt_s.*time_end_s./sum(dt_s);
        t_sill_count    = 0;
        dt_vr_count     = 0;

        while time<dt_vr % Time loop for heat conduction after sill emplacement
            if time_end_s-time>1e-6
                t_sill_count    = t_sill_count + 1;
                dt              = dt_sill(t_sill_count)*yr;            
            elseif dt_vr-time>1e-6
                dt_vr_count     = dt_vr_count + 1;
                dt              = dt_sill(end)*exp(dt_vr_count/3)*yr;
            end
            
            if time + dt/yr >= dt_vr
                dt = yr*(dt_vr - time);
            end
            
            time        = time + dt/yr;
            time_all    = time_all + dt/yr;

            
            %----------------------------------------------------------------------
            % Recompute Densities and heat capacities for coefficient in thermal
            % solver
            Rho_solid               = (1-Phi_rock).*Rho_rock + Phi_rock.*FLUID.rho;
            Cp_solid                = (1-Phi_rock).*Cp_rock + Phi_rock.*FLUID.cp;
            
            % Change heat capacity if sill T is between liquidus and solidus
            solid_ind               = find(Ind<0 & Temp<=Sol_rock);
            Cp_solid(solid_ind)     = Cp_rock(solid_ind);
            melt_solid              = find(Ind<0 & Temp>Sol_rock & Temp<Liq_rock);
            Rho_solid(melt_solid)   = Rho_melt(melt_solid);
            Cp_solid(melt_solid)    = Cp_melt(melt_solid).*(1 + (Lc_rock(melt_solid)./((Liq_rock(melt_solid)-Sol_rock(melt_solid)).*Cp_melt(melt_solid))));
            melt_ind                = find(Ind<0 & Temp>=Liq_rock);
            Rho_solid(melt_ind)     = Rho_melt(melt_ind);
            Cp_solid(melt_ind)      = Cp_melt(melt_ind);
            
            Coeff_temp              = Rho_solid.*Cp_solid;
            
            
            Temp_old = Temp;
            
            % -         Solve for Diffusion
            T = thermal1d_fem(Temp, Gcoord_n, Coeff_temp, K_solid, T_top, T_bot, dt, n_nel, n_nel+1, Latent_dehyd, Latent_om,...
                L2G, E2N);
            
            % Reorder Temperature in Global Mesh
            T_check(E2N(:)) = T(L2G(:));
            Temp(E2N(:)) = T(L2G(:));
            
            T_max(Temp>T_max) = Temp(Temp>T_max);
            
            %Calculate CO2 release from carbonates
            if max(ROCK.Carb)>0
                C1_ind      = find(Ind_carb(Active_nodes)==1);
                C2_ind      = find(Ind_carb(Active_nodes)==2);
                C3_ind      = find(Ind_carb(Active_nodes)==3);
                
                if ~isempty(C1_ind)
                    C1_now  = (1-Phi_rock(Active_nodes(C1_ind))).*Rho_rock(Active_nodes(C1_ind)).*(C1(Temp(Active_nodes(C1_ind)), Pres(Active_nodes(C1_ind)))./100);
                    C1_do  = find(CO2_release(Active_nodes(C1_ind))<C1_now);
                    CO2_release(Active_nodes(C1_ind(C1_do))) = C1_now(C1_do);
                end
                if ~isempty(C2_ind)
                    C2_now = (1-Phi_rock(Active_nodes(C2_ind))).*Rho_rock(Active_nodes(C2_ind)).*(C2(Temp(Active_nodes(C2_ind)), Pres(Active_nodes(C2_ind)))./100);
                    C2_do  = find(CO2_release(Active_nodes(C2_ind))<C2_now);
                    CO2_release(Active_nodes(C2_ind(C2_do))) = C2_now(C2_do);
                end
                if ~isempty(C3_ind)
                    C3_now = (1-Phi_rock(Active_nodes(C3_ind))).*Rho_rock(Active_nodes(C3_ind)).*(C3(Temp(Active_nodes(C3_ind)), Pres(Active_nodes(C3_ind)))./100);
                    C3_do  = find(CO2_release(Active_nodes(C3_ind))<C3_now);
                    CO2_release(Active_nodes(C3_ind(C3_do))) = C3_now(C3_do);
                end
            end
            
            % -  Calculate maturation
            [Latent_dehyd(Active_nodes), Latent_om(Active_nodes), Toc_rock_o(Active_nodes), W(Active_nodes,:), Ro(Active_nodes)] = coeff_diff(Temp(Active_nodes),...
                dt, Rho_rock(Active_nodes), A, R, length(Active_nodes), Ind(Active_nodes), Toc_rock(Active_nodes),Lm_rock(Active_nodes), Temp_old(Active_nodes),...
                Toc_rock_o(Active_nodes), Phi_rock(Active_nodes), W(Active_nodes,:), Ld_rock(Active_nodes));
            
            CO2_org(Active_nodes)     = (Toc_rock(Active_nodes)- Toc_rock_o(Active_nodes)).*(1-Phi_rock(Active_nodes)).*Rho_rock(Active_nodes).*3.66;
            
            Ro_plot = Ro;
            Ro_plot(Ind<0) = NaN;
            
            A_store                 = zeros(nnod,1);
            A_store(Active_nodes)   = 1;
            
            % Store Dynamic Variables
            Results.Temp(:,c1)            = Temp;                % Temperature
            Results.Pres(:,c1)            = Pres;                % Pressure
            Results.Toc(:,c1)             = Toc_rock_o.*1e2;     % TOC contents
            Results.Ro(:,c1)              = Ro_plot;             % Maturity
            Results.Tmax(:,c1)            = T_max;               % Maximum Temperature
            Results.Active(:,c1)          = A_store;             % Active nodes for time step
            Results.Gcoord(:,c1)          = Gcoord_n;            % Z-values for lithologies
            Results.CO2_release(:,c1)     = CO2_release;         % CO2 release
            Results.Time(c1,1)            = time_all;            % Time
            Results.CO2_org(:,c1)         = CO2_org;             % TOC reacted to CO2
            
            c1 = c1 + 1;
        end
        
    elseif ~isempty(ero)
        
        %% - Erosion
        % Remove layer that is eroded
        Ind_ero                     = find(ROCK.Ero_t==Times_all(i));
        for q = 1:length(Ind_ero)
            top     = ROCK.Tops(Ind_ero(q));
            if Ind_ero(q)<max(Ind_nel)
                bot     = ROCK.Tops(Ind_ero(q))+ROCK.Ero_thick(Ind_ero(q));
            else
                bot     = max(Gcoord_n);
            end
            ss = min(Gcoord_n(Ind==Ind_ero(q)));
            Gcoord_n(Gcoord_n>ss & Gcoord_n<=ss+ROCK.Ero_thick(Ind_ero(q))) = ss;
            Gcoord_n(Gcoord_n>ss) = Gcoord_n(Gcoord_n>ss) - (bot-top);
            Active(Gcoord_c>top & Gcoord_c<bot) = 0;
        end
        
        %         Set time allowed for VR calculation
        if i < length(Times_all)
            dt_vr         = (Times_all(i)-Times_all(i+1))*1e6*yr;
        else
            dt_vr         = 1e-2;
        end
        
        time_all = time_all + dt_vr/yr;
        
        % Create new mesh only for thermal diffusion
        n_nel = size(ELEM2NODE(:,Active'),2);
        E2N   = ELEM2NODE(:,Active');
        L2G = [1:n_nel; 2:n_nel+1];
        Active_nodes = unique(E2N(:));
        
        Pres(Active_nodes(1))       = min(Gcoord).*1.05e3*9.81;
        Pres(Active_nodes(2:end))   = cumsum(Rho_rock(Active_nodes(1:end-1)).*diff(Gcoord_n(Active_nodes)).*9.81) + Pres(Active_nodes(1));
        
        T_top                   = WellData.T(1,2);
        T_bot                   = (max(Gcoord_n(E2N(:)))).*polya(1)+polya(2);
        
        % Recompute Densities and heat capacities for coefficient in thermal
        % solver
        Rho_solid               = (1-Phi_rock).*Rho_rock + Phi_rock.*FLUID.rho;
        Cp_solid                = (1-Phi_rock).*Cp_rock + Phi_rock.*FLUID.cp;
        
        % Change heat capacity if sill T is between liquidus and solidus
        solid_ind               = find(Ind<0 & Temp<=Sol_rock);
        Cp_solid(solid_ind)     = Cp_rock(solid_ind);
        melt_solid              = find(Ind<0 & Temp>Sol_rock & Temp<Liq_rock);
        Rho_solid(melt_solid)   = Rho_melt(melt_solid);
        Cp_solid(melt_solid)    = Cp_melt(melt_solid).*(1 + (Lc_rock(melt_solid)./((Liq_rock(melt_solid)-Sol_rock(melt_solid)).*Cp_melt(melt_solid))));
        melt_ind                = find(Ind<0 & Temp>=Liq_rock);
        Rho_solid(melt_ind)     = Rho_melt(melt_ind);
        Cp_solid(melt_ind)      = Cp_melt(melt_ind);
        
        Coeff_temp              = Rho_solid.*Cp_solid;
        
        
        Temp_old = Temp;
        T_max(Temp>T_max) = Temp(Temp>T_max);
        
        % Solve for Diffusion
        T = thermal1d_fem(Temp, Gcoord_n, Coeff_temp, K_solid, T_top, T_bot, dt_vr, n_nel, n_nel+1, Latent_dehyd, Latent_om,...
            L2G, E2N);
        
        % Reorder Temperature in Global Mesh
        T_check(E2N(:)) = T(L2G(:));
        Temp(E2N(:))    = T(L2G(:));
        
        T_max(Temp>T_max) = Temp(Temp>T_max);
        
        %Calculate CO2 release from carbonates
        if max(ROCK.Carb)>0
            C1_ind      = find(Ind_carb(Active_nodes)==1);
            C2_ind      = find(Ind_carb(Active_nodes)==2);
            C3_ind      = find(Ind_carb(Active_nodes)==3);
            
            if ~isempty(C1_ind)
                C1_now  = (1-Phi_rock(Active_nodes(C1_ind))).*Rho_rock(Active_nodes(C1_ind)).*(C1(Temp(Active_nodes(C1_ind)), Pres(Active_nodes(C1_ind)))./100);
                C1_do  = find(CO2_release(Active_nodes(C1_ind))<C1_now);
                CO2_release(Active_nodes(C1_ind(C1_do))) = C1_now(C1_do);
            end
            if ~isempty(C2_ind)
                C2_now = (1-Phi_rock(Active_nodes(C2_ind))).*Rho_rock(Active_nodes(C2_ind)).*(C2(Temp(Active_nodes(C2_ind)), Pres(Active_nodes(C2_ind)))./100);
                C2_do  = find(CO2_release(Active_nodes(C2_ind))<C2_now);
                CO2_release(Active_nodes(C2_ind(C2_do))) = C2_now(C2_do);
            end
            if ~isempty(C3_ind)
                C3_now = (1-Phi_rock(Active_nodes(C3_ind))).*Rho_rock(Active_nodes(C3_ind)).*(C3(Temp(Active_nodes(C3_ind)), Pres(Active_nodes(C3_ind)))./100);
                C3_do  = find(CO2_release(Active_nodes(C3_ind))<C3_now);
                CO2_release(Active_nodes(C3_ind(C3_do))) = C3_now(C3_do);
            end
        end
        
        % Calculate maturation
        [Latent_dehyd(Active_nodes), Latent_om(Active_nodes), Toc_rock_o(Active_nodes), W(Active_nodes,:), Ro(Active_nodes)] = coeff_diff(Temp(Active_nodes),...
            dt_vr, Rho_rock(Active_nodes), A, R, length(Active_nodes), Ind(Active_nodes), Toc_rock(Active_nodes),Lm_rock(Active_nodes), Temp_old(Active_nodes),...
            Toc_rock_o(Active_nodes), Phi_rock(Active_nodes), W(Active_nodes,:), Ld_rock(Active_nodes));
        
        CO2_org(Active_nodes)     = (Toc_rock(Active_nodes)- Toc_rock_o(Active_nodes)).*(1-Phi_rock(Active_nodes)).*Rho_rock(Active_nodes).*3.66;
        
        Ro_plot = Ro;
        Ro_plot(Ind<0) = NaN;
        
        A_store                 = zeros(nnod,1);
        A_store(Active_nodes)   = 1;
        
        % Store Dynamic Variables
        Results.Temp(:,c1)            = Temp;                % Temperature
        Results.Pres(:,c1)            = Pres;                % Pressure
        Results.Toc(:,c1)             = Toc_rock_o.*1e2;     % TOC contents
        Results.Ro(:,c1)              = Ro_plot;             % Maturity
        Results.Tmax(:,c1)            = T_max;               % Maximum Temperature
        Results.Active(:,c1)          = A_store;             % Active nodes for time step
        Results.Gcoord(:,c1)          = Gcoord_n;            % Z-values for lithologies
        Results.CO2_release(:,c1)     = CO2_release;         % CO2 release
        Results.Time(c1,1)            = time_all;            % Time
        Results.CO2_org(:,c1)         = CO2_org;             % TOC reacted to CO2
        
        c1 = c1 + 1;
    end
    
end % End All Events
textprogressbar(' Done');

%% Compute Time Integrated Variables
% Use full column height to include eroded layers
Full_col     = Results.Gcoord(:,1);
Full_col_mat = repmat(Full_col, 1, c1-1);
Time.CO2_org = (Results.CO2_org(1:end-1,:) + Results.CO2_org(2:end,:))./2.*diff(Full_col_mat);
Time.CO2_rel = (Results.CO2_release(1:end-1,:) + Results.CO2_release(2:end,:))./2.*diff(Full_col_mat);