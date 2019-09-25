function [] = H_atom_DC_Stark_resonance
%-----------------------------------------------------------------
% Solve radial coupled equations for the DC Stark problem for H atom 
% Written by Tsogbayar Tsednee (PhD), Texas A&M University; 
% Contact: tsog215@gmail.com
% Feb 14, 2017
% ----------------------------------------------------------------
clear;
clc;
format long
N = 256; a = 0.; b = 100.; % Initial computaional parameters; you may change them
[r,wr,D]=legDC2(N,a,b);   % Legendre diff. matrix D; radial coordinate r; weight wr
D1 = (2/(b-a))*D; D2 = (2/(b-a))^2*D^2; % differentation matrix
Z = 1.; % charge
H0 = -0.5*D2 - diag((Z./r)); % S state Hamiltonian 
En = sort(eig(H0(2:N,2:N))); % Eigenvalue problem with boundary conditions
[En(1),En(2),En(3),En(4),En(5)] % eigenenergies for first 5 S-states
% Results
%  -0.500000000000000  -0.125000000000002  -0.055555555555557  -0.031249999999598  -0.019999971503745
%%%  --- Complex absorbing potential ---
eta = 0.021; % CAP strength parameter
ci = sqrt(-1.);
r_c = 70.00 ; % where CAP starts 
Vc = zeros(N,1);
for i = 1:N
    if (r(i) < r_c)
        Vc(i) = 0.; 
    else
        Vc(i) = -ci*eta*(r(i)-r_c).^2;
    end
end
Vc = Vc(2:N);
Vc = diag(Vc);
%%% for ell = 0, 1;
%%% DC Stark shift
F_str = 0.075; % DC field field strength in au; you may change it
H12 = F_str *diag(r(2:N)) * (1./sqrt(3.)); % off-diagonal dipole coupling elements
H21 = H12;
H11 = H0(2:N,2:N) + Vc;
%
ell = 1.
H_ham_ell_1 = -0.5*D2 + diag(ell*(ell+1)./(2.*r.^2)) - diag((Z./r));
H22 = H_ham_ell_1(2:N,2:N) + Vc;
H_ham_ell_1 = [H11, H12; 
               H21, H22];
%%% complex eigenvalue calculation begins
Ham_dc = sparse(H_ham_ell_1); % ell = 0, 1
rledc = -0.50; % 
opts.p = 48.;
%
En_res_par = eigs ( Ham_dc,5, rledc, opts ) % 
%%% complex eigenvalue calculation ends
%%%
ell = 12 % you may change it
%
delta_ell = eye(ell); % delta_ell matrix
%
delta_ellmellp1 = zeros(ell);
for i = 1:ell
    delta_ellmellp1(i,i) = (i-1)*((i-1)+1); % ell(ell+1) term
end
%
diag_term = zeros(ell); % diagonal elements from angular part 
for i = 2:ell-1
    diag_term(i,i+1) = sqrt(((i-1)+1)^2/((2*(i-1)+1)*(2*(i-1)+3)));
    diag_term(i,i-1) = sqrt(((i-2)+1)^2/((2*(i-2)+1)*(2*(i-2)+3)));     
end
diag_term(1,2) = sqrt(((1-1)+1)^2/((2*(1-1)+1)*(2*(1-1)+3)));
diag_term(ell,ell-1) = sqrt(((ell-2)+1)^2/((2*(ell-2)+1)*(2*(ell-2)+3)));
%
H_ham1 = kron(H0(2:N,2:N), delta_ell);  % 
H_ham2 = kron(diag(1./(2.*r(2:N).^2)), delta_ellmellp1);%  
H_ham3 = kron(Vc, delta_ell); % CAP term
H_ham4 = H_ham1 + H_ham2 + H_ham3; 
%
H_ham5 = F_str * kron(diag(r(2:N)), diag_term); % F*r*cos(theta)
%
H_ham_dc = sparse(H_ham4 + H_ham5);
%%% complex eigenvalue calculation begins
En_res_par = eigs (H_ham_dc,5, rledc, opts ) %  
%%% complex eigenvalue calculation ends  
% Results:
% -0.506105 -0.00003861i for F = 0.05 au
% -0.506105 -0.00003849i is from J. Phys. B: At. Mol. Opt. Phys. 33, 2195 (2000)
%%%
% -0.515257 -0.001503i for F = 0.075 au
% -0.515258 -0.001502i is from J. Phys. B: At. Mol. Opt. Phys. 33, 2195 (2000)
% %%
% -0.527416 -0.007271i for F = 0.1 au
% -0.527412 -0.007270i  is from J. Phys. B: At. Mol. Opt. Phys. 33, 2195 (2000)

