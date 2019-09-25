function [] = H2plus_eig_values_for_sigma_states
%-----------------------------------------------------------------
% Compute eigenvalue for the sigma states for H2+ ion in prolate
% spheroidal coordinate
% Written by Tsogbayar Tsednee (PhD), Texas A&M University; 
% Contact: tsog215@gmail.com
% Feb 14, 2017
% ----------------------------------------------------------------
clear;
clc;
N = 12; M = 10.; a = 1.; b = 20.; % Initial computaional parameters; you may change them
R = 2.0; % internuclear separation; you may change it
[mu,wrm,Dmu]=legDC2(N,a,b);  % Legendre diff. matrix and coordinate mu and weight
[nu,wrn,Dnu]=legDC2(M,-a,a); % Legendre diff. matrix and coordinate nu and weight
Dmu2 = Dmu*Dmu; Dnu2 = Dnu*Dnu; % kinetic energy matrix
Tmu = diag(mu.^2-1.)*(2/(b-a))^2*Dmu2 + 2.*diag(mu)*(2/(b-a))*Dmu;
Tnu = diag(1.-nu.^2)*Dnu2 - 2.*diag(nu)*Dnu;  
[mum,nun] = meshgrid(mu(2:N+1),nu); mum = mum(:); nun = nun(:);
Smunu= diag(mum.^2 - nun.^2);
Tmunu = (4./R^2)*(kron(Tmu(2:N+1,2:N+1),eye(M+1)) + ...
                     kron(eye(N),Tnu(1:M+1,1:M+1)));
Vcmunu = -(4./R)*diag(mum); % Coulomb potential 
Hmunu = -0.5*Tmunu + Vcmunu; % H0 hamiltonian matrix
  
En = sort(eig(Hmunu,Smunu));  % solve eigenvalues 
[En(1),En(2),En(3),En(4),En(5)] % eigenenergies for first 5 S-states
% Results
% -1.1024   -0.6676   -0.3609   -0.2554   -0.2358

