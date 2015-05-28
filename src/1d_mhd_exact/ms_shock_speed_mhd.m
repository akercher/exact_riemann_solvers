%******************************************************************************* 
%* Program: ms_shock_speed_mhd.m
%* Description: Computes speed of magnetosonic shock.  All variables are
%*              considered to be downstream (pre-shock) state except
%*              difference in total pressure.
%* Author: Andrew Kercher 
%* References: 
%*         [1] Dai, W. and Woodward, P. R., "An Approximate Riemann Solver
%*             for Ideal Magnetohydrodynamics", J. Comp. Phys.,
%*             111:354-373, 1994.
%-------------------------------------------------------------------------------
%******************************************************************************* 

% magnetosonic shock separating upstream (u) and downstream (d) regions
pt_d = pg_d + 0.5*(by_d^2 + bz_d^2);

% difference between upstream snd downstream total pressure
pt_ud = pt_u - pt_d;

% speed of sound
c0 = sqrt(gamma*pg_d/d_d);

% speeds in Lagrangian mass coordiantes
C0 = d_d*c0;
Cperp = sqrt(d_d*(by_d^2 + bz_d^2));

% Coefficients for magnetosonic shocks
Cms = sqrt(C0^2 + Cperp^2);
G1 = 0.5*(gamma + 3)*d_d*pt_ud;
G2 = 2*d_d^2*((gamma + 1)*pt_ud + 2*gamma*pt_d)*pt_ud;

% compute magnetosonic shocks speed
radical = sqrt((Cms^2 + G1)^2 - G2);
wms_2 = 0.5*(Cms^2 + G1 + radical);
wms = sqrt(wms_2);

