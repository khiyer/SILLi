function [Temp] = thermal1d_fem(Temp, Gcoord, coeff_temp, K_r, T_top, T_bot, dt, nel, nnod, Latent_dehyd,...
    Latent_om, L2G, NODES, t_bc, k_avg)
%  thermal1d_fem
%  
%  1D FEM code for thermal diffusion
%
%  Developed by Karthik Iyer, Henrik Svensen and Daniel W. Schmid
%
%--------------------------------------------------------------------------

%%
npe     = 2;
nip     = 2;

gauss   = [-sqrt(1/3) sqrt(1/3)];   
w       = [1 1];                                      %weights

KK_dif  = zeros(nel, npe*npe);
F_dif   = zeros(nnod,1);

KKi     = zeros(nel, npe*npe);
KKj     = zeros(nel, npe*npe);
for iel=1:nel
    cols        = ones(npe,1)*L2G(:,iel)';
    rows        = L2G(:,iel)*ones(1,npe);
    KKi(iel,:)  = rows(:);
    KKj(iel,:)  = cols(:);
end


for iel=1:nel
    AEL_dif = zeros(npe,npe);           %element stiffness matrix
    Fel_dif = zeros(npe,1);             %element force vector

    xeg     = Gcoord(NODES(:,iel));   %element x-coord. of global system

    Tel_dif = Temp(NODES(:,iel));
    E_d     = Latent_dehyd(NODES(:,iel));
    E_o     = Latent_om(NODES(:,iel));
    ED      = (0.5.*(K_r(NODES(1,iel)) + K_r(NODES(2,iel))));
    ER      = (0.5.*(coeff_temp(NODES(1,iel)) + coeff_temp(NODES(2,iel))));
    
    for ip=1:nip
        N           = [0.5*(1-gauss(ip)), 0.5*(1+gauss(ip))];

        dNdx        = [-1/(xeg(2) - xeg(1)), 1/(xeg(2) - xeg(1))];

        M_lumped    = diag(sum(N'*N,2));
        weight      = (xeg(2) - xeg(1))/2* w(ip);
        
        AEL_dif     = AEL_dif + (dt*ED*(dNdx' * dNdx) + M_lumped*ER) * weight;
        Fel_dif     = Fel_dif + (M_lumped*Tel_dif*ER - M_lumped*dt*E_d - M_lumped*dt*E_o)*weight;
    end


    % store element stiffness matrix
    KK_dif(iel,:) = AEL_dif(:);
    F_dif(L2G(:,iel))              = F_dif(L2G(:,iel))             + Fel_dif;
end

A_dif = sparse(KKi(:),KKj(:),KK_dif(:));

if t_bc.choice == 1
    % Apply bc (Constant T)
    Bc_it = [1, nnod];
    Bc_vt = [T_top, T_bot];
    
    for i=1:size(Bc_it,2)
        A_dif(Bc_it(i),:)        = 0;
        A_dif(Bc_it(i),Bc_it(i)) = 1;
        F_dif(Bc_it(i))          = Bc_vt(i);
    end
elseif t_bc.choice == 2
    % Apply bc (HF)
    Bc_it = 1;
    Bc_vt = T_top;
    
    for i=1:size(Bc_it,2)
        A_dif(Bc_it(i),:)        = 0;
        A_dif(Bc_it(i),Bc_it(i)) = 1;
        F_dif(Bc_it(i))          = Bc_vt(i);
    end
    F_dif(nnod)                  = F_dif(nnod) + t_bc.factor*k_avg*t_bc.t_grad*dt;
end

%Solution
Temp  =A_dif\F_dif;