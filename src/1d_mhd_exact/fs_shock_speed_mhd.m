%******************************************************************************* 
%* Program: fs_shock_speed_mhd.m
%* Description: Computes fast/slow shock speeds.  All variables are
%*              considered to be upstream (pre-shock) state except
%*              difference in bperp.
%* Author: Andrew Kercher 
%* References: 
%*         [1] Dai, W. and Woodward, P. R., "An Approximate Riemann Solver
%*             for Ideal Magnetohydrodynamics", J. Comp. Phys.,
%*             111:354-373, 1994.
%*         [2] Ryu, D. and Jones, T. W., "Numerical Magnetohydrodynamics in
%*             Astrophysics: Algorithm and Tests for One-Dimensional Flow",
%*             ApJ, 442:228-258, March 1995.
%-------------------------------------------------------------------------------
%******************************************************************************* 

% variable in mass coordinates
V_u = 1/d_u;

% Sound, Alfven, fast and slow characteristics speeds in Lagrangian mass coordiantes
[C02,Ca2,Cperp2,Cs2,Cf2] = Lagrangian_wave_speeds_mhd(gamma,d_u,pg_u,bperp_u,bx);

C0 = sqrt(C02);
Ca = sqrt(Ca2);
Cperp = sqrt(Cperp2);
Cf = sqrt(Cf2);
Cs = sqrt(Cs2);

if strcmp(method,'RJ')
%******************************************************************************* 
%* Using method of RJ [2]
%******************************************************************************* 

if min(bperp_u,bperp_d)/max(bperp_u,bperp_d) < tol
  if bperp_u < bperp_d
    bperp_u = tol*bperp_d;
    bperp_du = (1 - tol)*bperp_d;    
    bperp_r = bperp_du/(bperp_u);    
  else 
    bperp_d = tol*bperp_u;  
    bperp_du = (tol - 1)*bperp_u;    
    bperp_r = bperp_du/(bperp_u);    
  end
else
  % difference across shock of perpendicular magnetic field component  
  bperp_du = bperp_d - bperp_u;
  bperp_r = bperp_du/bperp_u;        
end  

% Coefficients, eq. (3.9) - (3.11) [2] (with correct def. of Cperp.)
S0 = -0.5*(gamma - 1)*bperp_r;

S1 = 0.5*(-(gamma-2)*Cperp2*bperp_r + 2*C02 ... 
	  - (gamma-4)*Cperp2 - 2*gamma*Ca2)*bperp_r;
     
S2 = 0.5*(Ca2*(bperp_du^2)/V_u + (gamma+2)*Cperp*Ca2*bperp_du/sqrt(V_u) ... 
	  + (gamma + 1)*Cperp2*Ca2 + (gamma + 1)*Ca2*Ca2 ... 
	  - 2*C02*Ca2)*bperp_r;     
     
% compute fast/slow shock speed     
[wfast2,wslow2,errio] = Lagrangian_shock_speeds_mhd(Cs2,Cf2,S0,S1,S2);
if errio == 1
  error('ShockSpeed:ImRadical','shock speed negative radical...');
end

wfast = sqrt(wfast2);
wslow = sqrt(wslow2);

end

if strcmp(method,'DW')
%******************************************************************************* 
%* Using method of DW [1]
%******************************************************************************* 
% speed in y-direction
Cy2 = (by_d^2)/V_d;
Cy = sqrt(Cy2);

% lambda defined and difference across shock of magnetic field y-component [1]
lambda_dw = 1 + (bz_u/by_u)^2;
by_du = by_d - by_u;
by_r = by_du/by_u;

% Coefficients from DW [1], S1 corrected to contain minus sign on first term
S0_dw = -0.5*(gamma - 1)*by_r;

S1_dw = 0.5*(-(gamma-2)*lambda_dw*by_u*by_du/V_u + 2*C02 ... 
	     - (gamma-4)*lambda_dw*Cy2 - 2*gamma*Ca2)*by_r;

S2_dw = 0.5*(lambda_dw*bx*bx*by_du*by_du/(V_u*V_u) ... 
	     + (gamma+2)*lambda_dw*bx*bx*by_u*by_du/(V_u*V_u) ... 
	     + lambda_dw*(gamma+1)*Ca2*Cy2 + (gamma + 1)*Ca2*Ca2 ... 
	     - 2*C02*Ca2)*by_r;     
	
% compute fast/slow shock speed     
[wfast_dw,wslow_dw] = Lagrangian_shock_speeds_mhd(Cs2,Cf2,S0_dw,S1_dw,S2_dw);

wfast = wfast_dw;
wslow = wslow_dw;
%* fprintf('RJ: wfast = %f wslow = %f\n',wfast,wslow);
%* fprintf('DW: wfast = %f wslow = %f\n',wfast_dw,wslow_dw);

end
