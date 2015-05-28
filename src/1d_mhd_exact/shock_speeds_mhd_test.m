%******************************************************************************* 
%* Program: shock_speeds_mhd_test.m
%* Description: Computes characteristic speeds for plane symmetric MHD
%*              equations.  
%* Author: Andrew Kercher 
%* References: 
%*         [1] Dai, W. and Woodward, P. R., "An Approximate Riemann Solver
%*             for Ideal Magnetohydrodynamics", J. Comp. Phys.,
%*             111:354-373, 1994.
%*         [2] Ryu, D. and Jones, T. W., "Numerical Magnetohydrodynamics in
%*             Astrophysics: Algorithm and Tests for One-Dimensional Flow",
%*             ApJ, 442:228-258, March 1995.
%*         [3] Toro, E. F., "Riemann Solvers and Numerical Methods for
%*             Fluid Dynamics", Springer-Verlag, 2nd ed., 1999.
%-------------------------------------------------------------------------------
%******************************************************************************* 

%* clear all;

%* bx = 1.3;

%* d_u = 0.6820;
%* pg_u = 0.5290;
%* by_u = 0.05555;%70190;
%* bz_u = 0.0;
bperp_u = sqrt(by_u^2 + bz_u^2);

%* by_d = -0.3442;
%* bz_d = 0.0;
bperp_d = sqrt(by_d^2 + bz_d^2);

% variable in mass coordinates
V_u = 1/d_u;

% Sound, Alfven, fast and slow characteristics speeds in Lagrangian mass coordiantes
[C02,Ca2,Cperp2,Cs2,Cf2] = Lagrangian_wave_speeds_mhd(gamma,d_u,pg_u,bperp_u,bx);

C0 = sqrt(C02);
Ca = sqrt(Ca2);
Cperp = sqrt(Cperp2);
Cf = sqrt(Cf2);
Cs = sqrt(Cs2);

if (abs(bx) > 0)

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
  
  fprintf('RJ: wfast = %f wslow = %f\n',wfast,wslow);

  %mach number
  Maf = wfast/Cf;
  Mar = wfast/Ca;
  Mas = wslow/Cs;

  fprintf(' Maf = %f Mar = %f Mas = %f\n',Maf,Mar,Mas)

else

  % Coefficients for magnetosonic shocks
  Cms = sqrt(C0^2 + Cperp^2);
  pt_d = pt_up - pt;
  G1 = 0.5*(gamma + 3)*d*pt_d;
  G2 = 2*d^2*((gamma + 1)*pt_d + 2*gamma*pt)*pt_d;

  % compute magnetosonic shocks speed
  radical = sqrt((Cms^2 + G1)^2 - G2);
  wms_2 = 0.5*(Cms^2 + G1 + radical);
  wms = sqrt(wms_2);
  wmsr = wms;  
  wmsl = -wms;

  fprintf('DW: wmsl = %f wmsr = %f\n',wmsl,wmsr);

  %mach number
  Ma = wmsr/Cms;
  
  fprintf(' Ma = %f\n',Ma)  

end