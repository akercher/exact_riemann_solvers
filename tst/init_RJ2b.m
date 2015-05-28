%******************************************************************************* 
%* Program: init_RJ2b.m
%* Description: Inititialization file for problem 2b [1].
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
prob_descrip = 'Riemann problem from Figure 2b of Ryu & Jones (1995).';
author  = 'D. Ryu & T.W. Jones';
journal = 'ApJ 442, 228-258 (1995)';

%-------------------------------------------------------------------------------
% Change default values
%-------------------------------------------------------------------------------
delx = 1.0e-1;
tol = 1.0e-12;

%-------------------------------------------------------------------------------
% Define problem
%-------------------------------------------------------------------------------
%* run_type = 'cont';
state_dimensions = 'none';       % ICs nondimensionalized?
method = 'RJ';				% RJ or DW
derivative_method = 'RJ';		% RJ or DW
Jacobian_method = 'approx';		% nested or approx
init_guess = 'RJ';			% user or RJ

%-------------------------------------------------------------------------------
% plotting
%-------------------------------------------------------------------------------

% y-axis limits for plotting
ylim_d = [-0.2 3.2];  
ylim_vx = [-2 0.2];  
ylim_vy = [-1.7 2.3];  
ylim_vz = [-0.1 1.2];  
ylim_pg = [0 11];  
ylim_en = [1 21];  
ylim_by = [0 3.5];  
ylim_bz = [-0.05 .55];  
ylim_psi = [-2 14];  

%* plot1 = 'plot_mhd';

%-------------------------------------------------------------------------------
% paramaters
%-------------------------------------------------------------------------------
tf = 0.035;
gamma = 5/3;				% ratio of specific heats
bx = 2/sqrt(4*pi);
  
%-------------------------------------------------------------------------------
% initial values
%-------------------------------------------------------------------------------
% initial left state  
dl = 1;
vxl = 0;
vyl = 0;
vzl = 0;
pgl = 1;
byl = 6/sqrt(4*pi);
bzl = 0;
kel = 0.5*(vxl*vxl + vyl*vyl + vzl*vzl);
pbl = 0.5*(bx*bx + byl*byl + bzl*bzl);
enl = pgl/(gamma - 1) + dl*kel + pbl;
ptl = pgl + pbl;

% initial right state  
dr = 0.1;
vxr = 0;
vyr = 2;
vzr = 1;
pgr = 10;
byr = 1/sqrt(4*pi);
bzr = 0;
ker = 0.5*(vxr*vxr + vyr*vyr + vzr*vzr);
pbr = 0.5*(bx*bx + byr*byr + bzr*bzr);
enr = pgr/(gamma - 1) + dr*ker + pbr;
ptr = pgr + pbr;

% guess bperp in regions 2,4,7 and rotation angle in region 3.
by2 = 3.11018112046060;
bz2 = 0;  
  
by3 = 3.04845572425575;
bz3 = 0.616558431393090;  
  
by4 = 1.21054075568080;
bz4 = 0.244835148341273;

by7 = 0.247876982919187;
bz7 = 0;
  
by2 = 2.5;
bz2 = 0;  
  
by3 = 2.0;
bz3 = 5.0e-1;  
  
by4 = 1.0;
bz4 = 3.0e-1;

by7 = 3.5e-1;
bz7 = 0;

bp2 = sqrt(by2*by2 + bz2*bz2);  
bp4 = sqrt(by4*by4 + bz4*bz4);    
bp7 = sqrt(by7*by7 + bz7*bz7);
  
% guess rotation angle inbetween left and right rotational discontinuity
psi3 = atan2(bz3,by3);
