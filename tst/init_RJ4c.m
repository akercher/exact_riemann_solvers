%******************************************************************************* 
%* Program: init_RJ4c.m
%* Description: Inititialization file for problem 4c [1].
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
prob_descrip = 'Riemann problem from Figure 4c of Ryu & Jones (1995).';
author  = 'D. Ryu & T.W. Jones';
journal = 'ApJ 442, 228-258 (1995)';

%-------------------------------------------------------------------------------
% Change default values
%-------------------------------------------------------------------------------
delx = 1.0e-8;
tol = 1.0e-10;

%-------------------------------------------------------------------------------
% Define problem
%-------------------------------------------------------------------------------
%* run_type = 'cont';
state_dimensions = 'dimensional';       % ICs nondimensionalized?
method = 'RJ';				% RJ or DW
derivative_method = 'RJ';		% RJ or DW
Jacobian_method = 'approx';		% nested or approx
init_guess = 'hlld';			% user or RJ

%-------------------------------------------------------------------------------
% plotting
%-------------------------------------------------------------------------------

% y-axis limits for plotting
ylim_d = [0.5 1.1];  
ylim_vx = [0.35 0.7];  
ylim_vy = [-1.1 -0.21];  
ylim_vz = [-0.5 0.5];  
ylim_pg = [0.45 0.85];  
ylim_en = [1.3 2.2];  
ylim_by = [-0.1 0.65];  
ylim_bz = [-0.5 0.5];  
ylim_psi = [-100 100];  

%* plot1 = 'plot_mhd';

%-------------------------------------------------------------------------------
% paramaters
%-------------------------------------------------------------------------------
tf = 0.15;
gamma = 5/3;				% ratio of specific heats
bx = 0.75;
  
%-------------------------------------------------------------------------------
% initial values
%-------------------------------------------------------------------------------
% initial left state  
dl = 0.65;
vxl = 0.667;
vyl = -0.257;
vzl = 0;
pgl = 0.5;
byl = 0.55;
bzl = 0;
kel = 0.5*(vxl*vxl + vyl*vyl + vzl*vzl);
pbl = 0.5*(bx*bx + byl*byl + bzl*bzl);
enl = pgl/(gamma - 1) + dl*kel + pbl;
ptl = pgl + pbl;

% initial right state  
dr = 1.0;
vxr = 0.4;
vyr = -0.94;
vzr = 0;
pgr = 0.75;
byr = 1.0e-5;
bzr = 0;
ker = 0.5*(vxr*vxr + vyr*vyr + vzr*vzr);
pbr = 0.5*(bx*bx + byr*byr + bzr*bzr);
enr = pgr/(gamma - 1) + dr*ker + pbr;
ptr = pgr + pbr;

% guess bperp in regions 2,4,7 and rotation angle in region 3.
by2 = 5.6490e-1;
bz2 = 0.0;    
  
by3 = 5.6490e-1;
bz3 = 0.0;
  
by4 = 1.8817e-4;
bz4 = 0.0;  

by7 = 1.066e-5;
bz7 = 0.0;

bp2 = sqrt(by2*by2 + bz2*bz2);  
bp4 = sqrt(by4*by4 + bz4*bz4);    
bp7 = sqrt(by7*by7 + bz7*bz7);
  
% guess rotation angle inbetween left and right rotational discontinuity
psi3 = atan2(bz3,by3);
