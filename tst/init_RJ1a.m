%******************************************************************************* 
%* Program: init_RJ1a.m
%* Description: Inititialization file for problem 1a [1].
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
prob_descrip = 'Riemann problem from Figure 1a of Ryu & Jones (1995).';
author  = 'D. Ryu & T.W. Jones';
journal = 'ApJ 442, 228-258 (1995)';

%-------------------------------------------------------------------------------
% Change default values
%-------------------------------------------------------------------------------
tol = 1.0e-12;

%-------------------------------------------------------------------------------
% Define problem
%-------------------------------------------------------------------------------
%* run_type = 'cont';
state_dimensions = 'none';       % dimensional or none
method = 'RJ';				% RJ or DW
derivative_method = 'RJ';		% RJ or DW
Jacobian_method = 'approx';		% nested or approx
init_guess = 'hlld';			% user or RJ

%-------------------------------------------------------------------------------
% plotting
%-------------------------------------------------------------------------------

% y-axis limits for plotting
ylim_d = [0.5 4.5];
ylim_vx = [-12 12];  
ylim_vy = [-0.5 0.5];  
ylim_vz = [-0.5 0.5];  
ylim_pg = [-10 160];  
ylim_en = [40 250];  
ylim_by = [0 7];  
ylim_bz = [-2 2];  
ylim_psi = [-100 100];  

%-------------------------------------------------------------------------------
% paramaters
%-------------------------------------------------------------------------------
tf = 0.08;
gamma = 5/3;				% ratio of specific heats
bx = 5/sqrt(4*pi);
bperpl = 5/sqrt(4*pi);
bperpr = 5/sqrt(4*pi);
alphal = 0;
alphar = 0;

%-------------------------------------------------------------------------------
% initial values
%-------------------------------------------------------------------------------
% initial left state  
dl = 1.0;
vxl = 10;
vyl = 0;
vzl = 0;
pgl = 20;
byl = bperpl*cos(alphal);
bzl = bperpl*sin(alphal);
kel = 0.5*(vxl*vxl + vyl*vyl + vzl*vzl);
pbl = 0.5*(bx*bx + byl*byl + bzl*bzl);
enl = pgl/(gamma - 1) + dl*kel + pbl;
ptl = pgl + pbl;

% initial right state  
dr = 1.0;
vxr = -10;
vyr = 0;
vzr = 0;
pgr = 1;
byr = bperpr*cos(alphar);
bzr = bperpr*sin(alphar);
ker = 0.5*(vxr*vxr + vyr*vyr + vzr*vzr);
pbr = 0.5*(bx*bx + byr*byr + bzr*bzr);
enr = pgr/(gamma - 1) + dr*ker + pbr;
ptr = pgr + pbr;

% guess bperp in regions 2,4,7 and rotation angle in region 3.
by2 = 3.8389;
bz2 = 0;  
  
by3 = 3.8389;
bz3 = 0;  
  
by4 = 4.0380;
bz4 = 0;

by7 = 5.4272;
bz7 = 0;
  
by2 = 3.83883482140749;
bz2 = 0;  
  
by3 = 3.83883482140749;
bz3 = 0;  
  
by4 = 4.03786851239782;
bz4 = 0;

by7 = 5.42712401583788;
bz7 = 0;
  
by2 = 4.0;
bz2 = 0;  
  
by3 = 4.0;
bz3 = 0;  
  
by4 = 4.5;
bz4 = 0;

by7 = 5.0;
bz7 = 0;  
  
bp2 = sqrt(by2*by2 + bz2*bz2);  
bp4 = sqrt(by4*by4 + bz4*bz4);    
bp7 = sqrt(by7*by7 + bz7*bz7);
  
% guess rotation angle inbetween left and right rotational discontinuity
psi3 = atan2(bz3,by3);
