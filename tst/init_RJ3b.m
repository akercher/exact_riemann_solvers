%******************************************************************************* 
%* Program: init_RJ3b.m
%* Description: Inititialization file for problem 3b [1].
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
prob_descrip = 'Riemann problem from Figure 5a of Ryu & Jones (1995).';
author  = 'D. Ryu & T.W. Jones';
journal = 'ApJ 442, 228-258 (1995)';

%-------------------------------------------------------------------------------
% Define problem
%-------------------------------------------------------------------------------
state_dimensions = 'dimensional';       % ICs nondimensionalized?
sampling = 'on';       			% sample data?
method = 'RJ';				% RJ or DW
derivative_method = 'RJ';		% RJ or DW
Jacobian_method = 'approx';		% nested or approx
init_guess = 'RJ';			% user or RJ

%-------------------------------------------------------------------------------
% plotting
%-------------------------------------------------------------------------------
% y-axis limits for plotting
ylim_d = [0.4 1.1];
ylim_vx = [-1.2 1.2];  
ylim_vy = [-1.2 1.2];  
ylim_vz = [-1.2 1.2];  
ylim_pg = [0.2 1.1];  
ylim_en = [0.2 2.8];  
ylim_by = [0.4 1.1];  
ylim_bz = [-0.5 0.5];    
ylim_psi = [-100 100];  

%-------------------------------------------------------------------------------
% paramaters
%-------------------------------------------------------------------------------
tf = 0.1;  
gamma = 5/3;				% ratio of specific heats
bx = 0;
  
%-------------------------------------------------------------------------------
% initial values
%-------------------------------------------------------------------------------
% initial left state  
dl = 1;
vxl = -1.0;
vyl = 0;
vzl = 0;
pgl = 1;
byl = 1.0;
bzl = 0;
kel = 0.5*(vxl*vxl + vyl*vyl + vzl*vzl);
pbl = 0.5*(bx*bx + byl*byl + bzl*bzl);
enl = pgl/(gamma - 1) + dl*kel + pbl;
ptl = pgl + pbl;  

% initial right state  
dr = 1;
vxr = 1;
vyr = 0;
vzr = 0;
pgr = 1;
byr = 1.0;
bzr = 0;
ker = 0.5*(vxr*vxr + vyr*vyr + vzr*vzr);
pbr = 0.5*(bx*bx + byr*byr + bzr*bzr);
enr = pgr/(gamma - 1) + dr*ker + pbr;
ptr = pgr + pbr;  

% guess bperp in regions 2 and 7.
by2 = 4.9638e-1;
bz2 = 0.0;

by4 = by2;
bz4 = bz2;

by7 = 4.9638e-1;
bz7 = 0.0;

bp2 = sqrt(by2*by2 + bz2*bz2);  
bp4 = bp2;  
bp7 = sqrt(by7*by7 + bz7*bz7);
