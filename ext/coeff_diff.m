function [Latent_dehyd, latent_om, TOC_all_new, W, Ro] = coeff_diff(Temp, dt, Rho_rock, A, R, nnod,...
    ind, TOC_ini, l_om, Temp_old, TOC_prev, Phi_all, W, l_d)
%   1D Thermal cooling of intrusive
%   COPYRIGHT(C) 2015, GeoModelling Solutions GmbH
%
%   VR CALCULATION
%
% Original author:    Karthik Iyer
% Last committed:     $Revision: 0 $
% Last changed by:    $Author: karthik $
% Last changed date:  $Date: 2012-03-07 16:05:55 +0100 (Wed, 07 Mar 2012) $
%--------------------------------------------------------------------------

%%
no_react = 20;

f_i = [0.03 0.03 0.04 0.04 0.05 0.05 0.06 0.04 0.04 0.07 0.06 0.06 0.06 0.05 0.05 0.04 0.03 0.02 0.02 0.01];
E_i = [142 151 159 167 176 184 192 201 209 218 226 234 243 251 259 268 276 285 293 301].*1e3;

T_AV = (Temp+Temp_old)./2 + 273.15;
F1   = zeros(nnod,1);
for r=1:no_react
    
    % Reaction Rate
    k       = A*exp(-E_i(r)/R./T_AV);
    
    % Update amount
    % If amount is used up then it is 0
    W(:,r)  = max(W(:,r).*exp(-k.*dt), 0);
    
    
    % Update Reaction Progress Recorder
    F1   = F1 + f_i(r).*(1-W(:,r));
end

% Evaluate synthetic vitrinite
Ro = exp(-1.6 + 3.7*F1);

TOC_all_new = TOC_ini - TOC_ini.*F1;

R_om = (1-Phi_all).*Rho_rock.*(TOC_prev - TOC_all_new);
latent_om = l_om.*R_om./dt;

Latent_dehyd = zeros(nnod,1);
T_ind = find(Temp>=350 & Temp<=650 & ind>0 & Temp>Temp_old);
if ~isempty(T_ind)
    Latent_dehyd(T_ind) = ((1-Phi_all(T_ind)).*Rho_rock(T_ind).*l_d(T_ind))./(650-350).*((Temp(T_ind) - Temp_old(T_ind))./dt);
end
Latent_dehyd(ind<0) = 0;