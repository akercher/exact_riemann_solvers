%******************************************************************************* 
%* Program: Lagrangian_shock_speed_derivative.m
%* Description: Computes derivatives of the various wave/shock speeds
%*              associated with the ideal MHD eqs.  The shock speed
%*              derivatives are needed to find the Jacobian matrix for
%*              the nonlinear MHD solver. 
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

% density
d_u = 1/V_u;

% Sound, Alfven, fast and slow characteristics speeds in Lagrangian mass coordiantes
[C02,Ca2,Cperp2,Cs2,Cf2] = Lagrangian_wave_speeds_mhd(gamma,d_u,pg_u,bperp_u,bx);

Cy2 = (by_u^2)/V_u;

radical = sqrt((C02 + Ca2 + Cperp2)^2 - 4*C02*Ca2);
C0 = sqrt(C02);
Ca = sqrt(Ca2);
Cperp = sqrt(Cperp2);
Cf = sqrt(Cf2);
Cs = sqrt(Cs2);

% calculate derivatives of speeds squared
[dC02,dCa2,dCperp2,dCs2,dCf2] = Lagrangian_wave_speed_derivative(gamma,V_u,pg_u,...
				bperp_u,bx,dV_u,dpg_u,dbperp_u,C02,Ca2,Cperp2);

dCy2 = (2*by_u/V_u)*dby_u - ((by_u*by_u)/(V_u*V_u))*dV_u;

% calculate derivative of Cperp from dCperp2 via chain rule
if abs(Cperp) > 0
  dCperp = dCperp2/(2*Cperp);  
else  
  dCperp = 0;
end

if strcmp(derivative_method,'RJ')  
%******************************************************************************* 
%* Using method of RJ [2]
%******************************************************************************* 
% perpendicular component of magnetic field
bperp_du = bperp_d - bperp_u;
bperp_r = bperp_du/bperp_u;

% check for switch-on shocks and adjust bperp to stabilize solution
if bperp_du <= 1e-12
  if abs(bperp_u) < abs(bperp_d)
    bperp_u = 1e-3*bperp_d;
  else
    bperp_d = 1e-3*bperp_u;  
  end
end  

%* if abs(bperp_du) > 0
% wave speed coefficients using method of RJ [2], eq. (3.9) - (3.11) [2]
% Coefficients, eq. (3.9) - (3.11) [2] (with correct def. of Cperp.)
S0 = -0.5*(gamma - 1)*bperp_r;

S1 = 0.5*(-(gamma-2)*Cperp2*bperp_r + 2*C02 ... 
	  - (gamma-4)*Cperp2 - 2*gamma*Ca2)*bperp_r;
     
S2 = 0.5*(Ca2*(bperp_du^2)/V_u + (gamma+2)*Cperp*Ca2*bperp_du/sqrt(V_u) ... 
	  + (gamma + 1)*Cperp2*Ca2 + (gamma + 1)*Ca2*Ca2 ... 
	  - 2*C02*Ca2)*bperp_r;     

% derivative of difference in upstream and downstean perpendicular magnetic field 
dbperp_du = dbperp_d - dbperp_u;
dbperp_r = (dbperp_d - (bperp_d/bperp_u)*dbperp_u)/bperp_u;
     
% derivatives of coefficients
dS0 = -0.5*(gamma - 1)*dbperp_r;

dS1 = 0.5*(-(gamma - 2)*(Cperp2*dbperp_r + bperp_r*dCperp2) ...
	+ 2*dC02 - (gamma - 4)*dCperp2 - 2*gamma*dCa2)*bperp_r ... 
        + (S1/bperp_r)*dbperp_r;

dS2 = 0.5*((2*Ca2*bperp_du/V_u)*dbperp_du + ((bperp_du^2)/(V_u))*dCa2 ... 
	- (Ca2*(bperp_du^2)/(V_u*V_u))*dV_u ... 
	+ (gamma + 2)*(Cperp*Ca2*dbperp_du/sqrt(V_u) ...
        + Cperp*bperp_du*dCa2/sqrt(V_u) ... 
	+ Ca2*bperp_du*dCperp/sqrt(V_u) ...
	- Cperp*Ca2*bperp_du*dV_u/(2*V_u^(3/2))) ...
	+ (gamma + 1)*(Cperp2*dCa2 + Ca2*dCperp2 + 2*Ca2*dCa2) ...
	- 2*(C02*dCa2 + Ca2*dC02))*bperp_r ...
	+ (S2/bperp_r)*dbperp_r;

% Calculate derivative of wave speed squared.
radical = sqrt((Cs2 + Cf2 + S1)^2 - 4*(1 + S0)*(Cs2*Cf2 - S2));
dradical = ((Cs2 + Cf2 + S1)^2)*(dCs2 + dCf2 + dS1) ...
    - 2*(1 + S0)*(Cs2*dCf2 + Cf2*dCs2 - dS2) ...
    - 2*dS0*(Cs2*Cf2 - S2);

dradical = dradical/radical;

dwf2 = 0.5*(1/(1 + S0))*(dCs2 + dCf2 + dS1 + dradical) ... 
       - (wfs*wfs/(1 + S0))*dS0;
dws2 = 0.5*(1/(1 + S0))*(dCs2 + dCf2 + dS1 - dradical) ... 
       - (wfs*wfs/(1 + S0))*dS0;
       
% use chain rule to calculate derivative       
dwf = dwf2/(2*wfs);
dws = dws2/(2*wfs);

%* else
%*   dwf = 0;
%*   dws = 0;
%* end

end

if strcmp(derivative_method,'DW')
%******************************************************************************* 
%* Using method of DW [1]
%******************************************************************************* 

% lambda defined and difference across shock of magnetic field y-component [1]
lambda_dw = 1 + (bz_d/by_d)^2;
by_du = by_u - by_d;
by_r = by_du/by_d;
%* lambda_dw2 = lambda_dw*lambda_dw;

%* if abs(by_du) > 0
% Coefficients from DW [1]
S0 = -0.5*(gamma - 1)*by_r;

S1 = 0.5*(-(gamma-2)*lambda_dw*by_d*by_du/V_u + 2*C02 ... 
	     - (gamma-4)*lambda_dw*Cy2 - 2*gamma*Ca2)*by_r;

S2 = 0.5*(lambda_dw*bx*bx*by_du*by_du/(V_u*V_u) ... 
	     + (gamma+2)*lambda_dw*bx*bx*by_d*by_du/(V_u*V_u) ... 
	     + lambda_dw*(gamma+1)*Ca2*Cy2 + (gamma + 1)*Ca2*Ca2 ... 
	     - 2*C02*Ca2)*by_r;     

%* S0 = -0.5*(gamma - 1)*by_r;
%* S1 = 0.5*(-(gamma-2)*d_d*lambda_dw*by_d*by_du + 2*C0^2 ... 
%* 	     - (gamma-4)*lambda_dw*Cy^2 - 2*gamma*Ca^2)*by_r;
%* S2 = 0.5*(lambda_dw*d_d^2*bx^2*by_du^2 ... 
%* 	     + (gamma+2)*lambda_dw*d_d^2*bx^2*by_d*by_du ... 
%* 	     + lambda_dw*(gamma+1)*Ca^2*Cy^2 + (gamma + 1)*Ca^4 ... 
%* 	     - 2*C0^2*Ca^2)*by_r;     

% derivative of differecne in upstream and downstean y-magnetic field 
dby_du = dby_d - dby_u;
dby_r = (dby_d - (by_d/by_u)*dby_u)/by_u;
dlambda = -(2*bz_u*bz_u/(by_u*by_u*by_u))*dby_u;
V_u3 = V_u*V_u*V_u;

% derivatives of coefficients
dS0 = -0.5*(gamma - 1)*dby_r;

dS1 = 0.5*((gamma - 2)*((lambda_dw*by_u/V_u)*dby_du ... 
	+ (lambda_dw*by_du/V_u)*dby_d + (by_u*by_du/V_u)*dlambda ... 
	- (lambda_dw*by_u*by_du/(V_u*V_u))*dV_u) + 2*dC02 ...
	- (gamma - 4)*(lambda_dw*dCy2 + Cy2*dlambda ...
	- 2*gamma*dCa2)  )*by_r + (S2/by_r)*dby_r;

dS2 = 0.5*(((2*lambda_dw*bx*bx*by_du/V_u)*dby_du ...
	+ (bx*by_du/(V_u*V_u))*dlambda - (lambda_dw*bx*by_du/V_u3)*dV_u) ...
	+ (gamma + 2)*((lambda_dw*bx*bx*by_u/(V_u*V_u))*dby_du ...
	+ (lambda_dw*bx*bx*by_du/(V_u*V_u))*dby_u ...
	+ (bx*bx*by_u*by_du/(V_u*V_u))*dlambda ...
	- (lambda_dw*bx*bx*by_u*by_du/(V_u3))*dV_u) ...
	+ (gamma + 1)*(lambda_dw*Ca2*dCy2 + lambda_dw*Cy2*dCa2 ... 
	+ Ca2*Cy2*dlambda) ...
	+ (gamma + 1)*(2*Ca2*dCa2) ...
	- 2*(C02*dCa2 + Ca2*dC02))*by_r + (S2/by_r)*dby_r;
	
		
% Calculate derivative of wave speed squared.
radical = sqrt((Cs2 + Cf2 + S1)^2 - 4*(1 + S0)*(Cs2*Cf2 - S2));
dradical = ((Cs2 + Cf2 + S1)^2)*(dCs2 + dCf2 + dS1) ...
    - 2*(1 + S0)*(Cs2*dCf2 + Cf2*dCs2 - dS2) ...
    - 2*dS0*(Cs2*Cf2 - S2);

dradical = dradical/radical;

dwf2_dw = 0.5*(1/(1 + S0))*(dCs2 + dCf2 + dS1 + dradical) ... 
       - (wfs*wfs/(1 + S0))*dS0;
dws2_dw = 0.5*(1/(1 + S0))*(dCs2 + dCf2 + dS1 - dradical) ... 
       - (wfs*wfs/(1 + S0))*dS0;
       
% use chain rule to calculate derivative       
dwf = dwf2_dw/(2*wfs);
dws = dws2_dw/(2*wfs);

%* else
%*   dwf = 0;
%*   dws = 0;
%* end

end