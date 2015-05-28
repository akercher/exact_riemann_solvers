%******************************************************************************* 
%* Program: init_RJ4b.m
%* Description: Inititialization file for problem 4b [1].
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
prob_descrip = 'Riemann problem from Figure 4b of Ryu & Jones (1995).';
author  = 'D. Ryu & T.W. Jones';
journal = 'ApJ 442, 228-258 (1995)';

%-------------------------------------------------------------------------------
% Change default values
%-------------------------------------------------------------------------------
delx = 1.0e-4;
tol = 1.0e-12;
maxit = 20;

%-------------------------------------------------------------------------------
% Define problem
%-------------------------------------------------------------------------------
%* state_dimensions = 'dimensional';       % ICs nondimensionalized?
method = 'RJ';				% RJ or DW
derivative_method = 'RJ';		% RJ or DW
Jacobian_method = 'approx';		% nested or approx
init_guess = 'user';			% user or RJ

%-------------------------------------------------------------------------------
% plotting
%-------------------------------------------------------------------------------

% y-axis limits for plotting
ylim_d = [0.3 1.1];  
ylim_vx = [-0.8 0.1];  
ylim_vy = [-0.1 1.1];  
ylim_vz = [-0.5 0.5];  
ylim_pg = [0.4 1.1];  
ylim_en = [1.8 3.0];  
ylim_by = [-0.1 1.1];  
ylim_bz = [-0.5 0.5];  
ylim_psi = [-10 190];  

%-------------------------------------------------------------------------------
% paramaters
%-------------------------------------------------------------------------------
tf = 0.15;
gamma = 5/3;				% ratio of specific heats
bx = 1.3;

%-------------------------------------------------------------------------------
% initial values
%-------------------------------------------------------------------------------
% initial left state  
dl = 0.4;
vxl = -0.66991;
vyl = 0.98263;
vzl = 0;
pgl = 0.52467;
byl = 0.0025293;
bzl = 0;
kel = 0.5*(vxl*vxl + vyl*vyl + vzl*vzl);
pbl = 0.5*(bx*bx + byl*byl + bzl*bzl);
enl = pgl/(gamma - 1) + dl*kel + pbl;
ptl = pgl + pbl;

% initial right state  
dr = 1.0;
vxr = 0;
vyr = 0;
vzr = 0;
pgr = 1.0;
byr = 1.0;
bzr = 0;
ker = 0.5*(vxr*vxr + vyr*vyr + vzr*vzr);
pbr = 0.5*(bx*bx + byr*byr + bzr*bzr);
enr = pgr/(gamma - 1) + dr*ker + pbr;
ptr = pgr + pbr;

% guess bperp in regions 2,4,7 and rotation angle in region 3.
by2 = 2.5305e-3;
bz2 = 0.0;    
  
by3 = 2.5305e-3;
bz3 = 0.0;
  
by4 = 2.5305e-3;
bz4 = 0.0;  

by7 = 2.5304e-3;
bz7 = 0.0;

bp2 = sqrt(by2*by2 + bz2*bz2);  
bp4 = sqrt(by4*by4 + bz4*bz4);    
bp7 = sqrt(by7*by7 + bz7*bz7);
  
% guess rotation angle inbetween left and right rotational discontinuity
psi3 = atan2(bz3,by3);
