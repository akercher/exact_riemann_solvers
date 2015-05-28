%******************************************************************************* 
%* Program: init_RJ4d.m
%* Description: Inititialization file for problem 4d [1].
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
prob_descrip = 'Riemann problem from Figure 4d of Ryu & Jones (1995).';
author  = 'D. Ryu & T.W. Jones';
journal = 'ApJ 442, 228-258 (1995)';

%-------------------------------------------------------------------------------
% Change default values
%-------------------------------------------------------------------------------
delx = 1.0e-6;
tol = 1.0e-12;
maxit = 50;

%-------------------------------------------------------------------------------
% Define problem
%-------------------------------------------------------------------------------
%* run_type = 'cont';
state_dimensions = 'dimensional';       % dimensional or none
method = 'RJ';				% RJ or DW
derivative_method = 'RJ';		% RJ or DW
Jacobian_method = 'approx';		% nested or approx
init_guess = 'hlld';			% user or RJ

%-------------------------------------------------------------------------------
% plotting
%-------------------------------------------------------------------------------

% y-axis limits for plotting
ylim_d = [0.1 1.2];  
ylim_vx = [-0.1 0.4];  
ylim_vy = [-0.1 0.9];  
ylim_vz = [-0.1 1.1];  
ylim_pg = [-0.1 1.2];  
ylim_en = [0.9 1.9];  
ylim_by = [-0.1 1.1];  
ylim_bz = [-0.1 0.5];  
ylim_psi = [-5 40];  

%-------------------------------------------------------------------------------
% paramaters
%-------------------------------------------------------------------------------
tf = 0.16;
gamma = 5/3;				% ratio of specific heats
bx = 0.70;
  
dl = 1.0;
vxl = 0;
vyl = 0;
vzl = 0;
pgl = 1.0;
byl = 1.0e-3;
bzl = 0;
kel = 0.5*(vxl*vxl + vyl*vyl + vzl*vzl);
pbl = 0.5*(bx*bx + byl*byl + bzl*bzl);
enl = pgl/(gamma - 1) + dl*kel + pbl;
ptl = pgl + pbl;

dr = 0.3;
vxr = 0;
vyr = 0;
vzr = 1;
pgr = 0.2;
byr = 1;
bzr = 0;
ker = 0.5*(vxr*vxr + vyr*vyr + vzr*vzr);
pbr = 0.5*(bx*bx + byr*byr + bzr*bzr);
enr = pgr/(gamma - 1) + dr*ker + pbr;
ptr = pgr + pbr;

%-------------------------------------------------------------------------------
% initial values
%-------------------------------------------------------------------------------
% initial left state  
dl = 1.0;
vxl = 0;
vyl = 0;
vzl = 0;
pgl = 1.0;
byl = 1.0e-3;
bzl = 0;
kel = 0.5*(vxl*vxl + vyl*vyl + vzl*vzl);
pbl = 0.5*(bx*bx + byl*byl + bzl*bzl);
enl = pgl/(gamma - 1) + dl*kel + pbl;
ptl = pgl + pbl;

% initial right state  
dr = 0.3;
vxr = 0;
vyr = 0;
vzr = 1;
pgr = 0.2;
byr = 1;
bzr = 0;
ker = 0.5*(vxr*vxr + vyr*vyr + vzr*vzr);
pbr = 0.5*(bx*bx + byr*byr + bzr*bzr);
enr = pgr/(gamma - 1) + dr*ker + pbr;
ptr = pgr + pbr;

% y-axis limits for plotting
ylim_d = [0.1 1.2];  
ylim_vx = [-0.1 0.4];  
ylim_vy = [-0.1 0.9];  
ylim_vz = [-0.1 1.1];  
ylim_pg = [-0.1 1.2];  
ylim_en = [0.9 1.9];  
ylim_by = [-0.1 1.1];  
ylim_bz = [-0.1 0.5];  
ylim_psi = [-5 40];  

% guess bperp in regions 2,4,7 and rotation angle in region 3.
by2 = 9.17758620157921e-04;
bz2 = 0.0;  
  
by3 = 8.04840212916876e-04;
bz3 = 4.41036411814364e-04;  

by4 = 0.658615974437010;  
bz4 = 0.360908440541998;    
  
by7 = 0.990679609370528;
bz7 = 0.0;  
  
by2 = 9.0e-4;
bz2 = 0.0;    
  
by3 = 8.0e-1;
bz3 = 5.0e-1;
  
by4 = 7.0e-1;
bz4 = 3.0e-1;  

by7 = 9.0e-1;
bz7 = 0.0;

bp2 = sqrt(by2*by2 + bz2*bz2);  
bp4 = sqrt(by4*by4 + bz4*bz4);    
bp7 = sqrt(by7*by7 + bz7*bz7);
  
% guess rotation angle inbetween left and right rotational discontinuity
psi3 = atan2(bz3,by3);
