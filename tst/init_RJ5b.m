%******************************************************************************* 
%* Program: init_RJ5b.m
%* Description: Inititialization file for problem 5b [1].
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
prob_descrip = 'Riemann problem from Figure 5b of Ryu & Jones (1995).';
author  = 'D. Ryu & T.W. Jones';
journal = 'ApJ 442, 228-258 (1995)';

%-------------------------------------------------------------------------------
% Change default values
%-------------------------------------------------------------------------------
tol = 1.0e-12;
nrp = 256;

%-------------------------------------------------------------------------------
% Define problem
%-------------------------------------------------------------------------------
state_dimensions = 'dimensional';       % ICs nondimensionalized?
method = 'RJ';				% RJ or DW
derivative_method = 'RJ';		% RJ or DW
Jacobian_method = 'approx';		% nested or approx
init_guess = 'hlld';			% user or RJ

%-------------------------------------------------------------------------------
% plotting
%-------------------------------------------------------------------------------

% y-axis limits for plotting
ylim_d = [0.2 1.1];  
ylim_vx = [-0.7 0.8];  
ylim_vy = [-1.7 0.2];  
ylim_vz = [-1.0 1.0];  
ylim_pg = [0.2 1.1];  
ylim_en = [1.2 3.3];  
ylim_by = [-1.2 1.2];  
ylim_bz = [-1.2 1.2];  
ylim_psi = [-10 190];  

%-------------------------------------------------------------------------------
% paramaters
%-------------------------------------------------------------------------------
tf = 0.16;
gamma = 5/3;				% ratio of specific heats
bx = 1.3;
alpha = pi;

%-------------------------------------------------------------------------------
% initial values
%-------------------------------------------------------------------------------
% initial left state  
dl = 1.0;
vxl = 0;
vyl = 0;
vzl = 0;
pgl = 1.0;
byl = 1.0;%0.031150;%1.0;
bzl = 0.0;%0.34672;%0;
kel = 0.5*(vxl*vxl + vyl*vyl + vzl*vzl);
pbl = 0.5*(bx*bx + byl*byl + bzl*bzl);
enl = pgl/(gamma - 1) + dl*kel + pbl;
ptl = pgl + pbl;

% initial right state  
dr = 0.4;
vxr = 0;
vyr = 0;
vzr = 0;
pgr = 0.4;
byr = cos(alpha);
bzr = sin(alpha);
ker = 0.5*(vxr*vxr + vyr*vyr + vzr*vzr);
pbr = 0.5*(bx*bx + byr*byr + bzr*bzr);
enr = pgr/(gamma - 1) + dr*ker + pbr;
ptr = pgr + pbr;

% guess bperp in regions 2,4,7 and rotation angle in region 3.
by2 = 0.348504250010020;
bz2 = 0.0;  

by3 = -0.348504250010020;
bz3 = 0.0;

by4 = -0.237022654767933;  
bz4 = 0.0;
  
by7 = -0.441568383543591;  
bz7 = 0.0;  

by2 = 4.0e-1;
bz2 = 0.0;  
  
by3 = -3.0e-1;
bz3 = 0.0;  
  
by4 = -2.5e-1;
bz4 = 0.0;

by7 = -5.0e-1;
bz7 = 0.0;

bp2 = sqrt(by2*by2 + bz2*bz2);  
bp4 = sqrt(by4*by4 + bz4*bz4);    
bp7 = sqrt(by7*by7 + bz7*bz7);
  
% guess rotation angle inbetween left and right rotational discontinuity
psi3 = atan2(bz3,by3);

