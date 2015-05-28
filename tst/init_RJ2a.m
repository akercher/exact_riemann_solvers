%******************************************************************************* 
%* Program: init_RJ2a.m
%* Description: Inititialization file for problem 2a [1].
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
prob_descrip = 'Riemann problem from Figure 2a of Ryu & Jones (1995).';
author  = 'D. Ryu & T.W. Jones';
journal = 'ApJ 442, 228-258 (1995)';

%-------------------------------------------------------------------------------
% Change default values
%-------------------------------------------------------------------------------
tol = 1.0e-12;

%-------------------------------------------------------------------------------
% Define problem
%-------------------------------------------------------------------------------
state_dimensions = 'dimensional';       % ICs nondimensionalized?
state_dimensions = 'none';       % ICs nondimensionalized?
method = 'RJ';				% RJ or DW
derivative_method = 'RJ';		% RJ or DW
Jacobian_method = 'approx';		% nested or approx
init_guess = 'hlld';			% user or RJ

%-------------------------------------------------------------------------------
% plotting
%-------------------------------------------------------------------------------

% y-axis limits for plotting
ylim_d = [0.9 1.75];  
ylim_vx = [-0.1 1.3];  
ylim_vy = [-0.22 0.26];  
ylim_vz = [-0.1 0.6];  
ylim_pg = [0.8 2.1];  
ylim_en = [2 5];  
ylim_by = [0.9 1.7];  
ylim_bz = [0.35 0.85];  
ylim_psi = [15 31];  

%* plot1 = 'plot_mhd';

%-------------------------------------------------------------------------------
% paramaters
%-------------------------------------------------------------------------------
tf = 0.2;
gamma = 5/3;				% ratio of specific heats
bx = 2/sqrt(4*pi);

%-------------------------------------------------------------------------------
% initial values
%-------------------------------------------------------------------------------
% initial left state  
dl = 1.08;
vxl = 1.2;
vyl = 0.01;
vzl = 0.5;
pgl = 0.95;
byl = 3.6/sqrt(4*pi);
bzl = 2/sqrt(4*pi);
kel = 0.5*(vxl*vxl + vyl*vyl + vzl*vzl);
pbl = 0.5*(bx*bx + byl*byl + bzl*bzl);
enl = pgl/(gamma - 1) + dl*kel + pbl;
ptl = pgl + pbl;

% initial right state  
dr = 1.0;
vxr = 0;
vyr = 0;
vzr = 0;
pgr = 1;
byr = 4/sqrt(4*pi);
bzr = 2/sqrt(4*pi);
ker = 0.5*(vxr*vxr + vyr*vyr + vzr*vzr);
pbr = 0.5*(bx*bx + byr*byr + bzr*bzr);
enr = pgr/(gamma - 1) + dr*ker + pbr;
ptr = pgr + pbr;

% guess bperp in regions 2,4,7 and rotation angle in region 3.
by2 = 1.43831583126146;
bz2 = 7.99064350700814e-1;  

by3 = 1.57164643419684;
bz3 = 0.487015131301542;

by4 = 1.41255895316497;
bz4 = 4.37717775469546e-1;

by7 = 1.50784400893008;
bz7 = 7.53922004465042e-1;
  
by2 = 1.0;
bz2 = 8.0e-1;    

by3 = 2.0;
bz3 = 1.0e-1;    
  
by4 = 1.5;
bz4 = 9.0e-1;

by7 = 1.0;
bz7 = 9.0e-1;  

bp2 = sqrt(by2*by2 + bz2*bz2);  
bp4 = sqrt(by4*by4 + bz4*bz4);    
bp7 = sqrt(by7*by7 + bz7*bz7);
  
% guess rotation angle inbetween left and right rotational discontinuity
psi3 = atan2(bz3,by3);

