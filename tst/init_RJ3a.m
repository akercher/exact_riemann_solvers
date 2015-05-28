%******************************************************************************* 
%* Program: init_RJ3a.m
%* Description: Inititialization file for problem 3a [1].
%* Author: Andrew Kercher 
%* References: 
%*         [1] Ryu, D. and Jones, T. W., "Numerical Magnetohydrodynamics in
%*             Astrophysics: Algorithm and Tests for One-Dimensional Flow",
%*             ApJ, 442:228-258, March 1995.
%-------------------------------------------------------------------------------
%******************************************************************************* 

%-------------------------------------------------------------------------------
% Problem Information
%-------------------------------------------------------------------------------
prob_descrip = 'Riemann problem from Figure 3a of Ryu & Jones (1995).';
author  = 'D. Ryu & T.W. Jones';
journal = 'ApJ 442, 228-258 (1995)';

%-------------------------------------------------------------------------------
% Change default values
%-------------------------------------------------------------------------------
delx = 1.0e-6;

%-------------------------------------------------------------------------------
% Define problem
%-------------------------------------------------------------------------------
%* run_type = 'cont';
state_dimensions = 'dimensional';       % ICs nondimensionalized?
sampling = 'on';       			% sample data?
method = 'RJ';				% RJ or DW
derivative_method = 'RJ';		% RJ or DW
Jacobian_method = 'approx';		% nested or approx
init_guess = 'user';			% user or RJ

%-------------------------------------------------------------------------------
% plotting
%-------------------------------------------------------------------------------
% y-axis limits for plotting
ylim_d = [-0.15 0.45];
ylim_vx = [-5 55];  
ylim_vy = [-5 5];  
ylim_vz = [-5 5];  
ylim_pg = [-10 90];  
ylim_en = [-20 265];  
ylim_by = [-1.5 1.5];  
ylim_bz = [-2.8 2.8];  
ylim_psi = [-130 80];  

%-------------------------------------------------------------------------------
% paramaters
%-------------------------------------------------------------------------------
tf = 0.01;  
gamma = 5/3;				% ratio of specific heats
bx = 0;

%-------------------------------------------------------------------------------
% initial values
%-------------------------------------------------------------------------------
% initial left state  
dl = 0.1;
vxl = 50.0;
vyl = 0;
vzl = 0;
pgl = 0.4;
byl = -1.0/sqrt(4*pi);
bzl = -2.0/sqrt(4*pi);
kel = 0.5*(vxl*vxl + vyl*vyl + vzl*vzl);
pbl = 0.5*(bx*bx + byl*byl + bzl*bzl);
enl = pgl/(gamma - 1) + dl*kel + pbl;
ptl = pgl + pbl;  

% initial right state  
dr = 0.1;
vxr = 0;
vyr = 0;
vzr = 0;
pgr = 0.2;
byr = 1.0/sqrt(4*pi);
bzr = 2.0/sqrt(4*pi);
ker = 0.5*(vxr*vxr + vyr*vyr + vzr*vzr);
pbr = 0.5*(bx*bx + byr*byr + bzr*bzr);
enr = pgr/(gamma - 1) + dr*ker + pbr;
ptr = pgr + pbr;  

% guess bperp in regions 2 and 7.
by2 = -1.0921;
bz2 = -2.1842;

by4 = by2;
bz4 = bz2;

by7 = 1.1014;
bz7 = 2.2029;

bp2 = sqrt(by2*by2 + bz2*bz2);  
bp4 = bp2;  
bp7 = sqrt(by7*by7 + bz7*bz7);
