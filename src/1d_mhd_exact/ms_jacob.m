%******************************************************************************* 
%* Program: ms_jacob.m
%* Description: Computes Jacobian (derivative matrix) for nonlinear MHD
%*              equations.  Used to find exact solution to 1D shocktube
%*              problems. 
%* Author: Andrew Kercher 
%* References: 
%*         [1] Dai, W. and Woodward, P. R., "An Approximate Riemann Solver
%*             for Ideal Magnetohydrodynamics", J. Comp. Phys.,
%*             111:354-373, 1994.
%*         [2] Ryu, D. and Jones, T. W., "Numerical Magnetohydrodynamics in
%*             Astrophysics: Algorithm and Tests for One-Dimensional Flow",
%*             ApJ, 442:228-258, March 1995.
%-------------------------------------------------------------------------------
%* 
%* Calculates the simplified version for magnetosonic shocks of jacobian
%* matrix for the system given by equations (11a) - (11d) in [1].
%* Following [2] the y-component of the magnetic field is replaced with
%* the perpendicular component.
%* 
%* Shock speeds: calculated in rstates_mhd.m
%*     wmsl = fast magnetosonic wave moving in negative direction (left).
%*     wmsr = fast magnetosonic wave moving in positive direction (right).
%*          
%******************************************************************************* 

% initialize jacobian
DJ = zeros(2,2);

% define state variables in regions 1,2,3,4.
V1 = 1/pstates(1,1);
vx1 = pstates(2,1);
vy1 = pstates(3,1);
vz1 = pstates(4,1);
pg1 = pstates(5,1);
by1 = pstates(6,1);
bz1 = pstates(7,1);
bp1 = bperp(1);

V2 = 1/pstates(1,2);
vx2 = pstates(2,2);
vy2 = pstates(3,2);
vz2 = pstates(4,2);
pg2 = pstates(5,2);
by2 = pstates(6,2);
bz2 = pstates(7,2);
bp2 = bperp(2);

%* V3 = 1/pstates(1,3);
%* vx3 = pstates(2,3);
%* vy3 = pstates(3,3);
%* vz3 = pstates(4,3);
%* pg3 = pstates(5,3);
%* by3 = pstates(6,3);
%* bz3 = pstates(7,3);
%* bp3 = bperp(3);
%* 
%* V4 = 1/pstates(1,4);
%* vx4 = pstates(2,4);
%* vy4 = pstates(3,4);
%* vz4 = pstates(4,4);
%* pg4 = pstates(5,4);
%* by4 = pstates(6,4);
%* bz4 = pstates(7,4);
%* bp4 = bperp(4);

V7 = 1/pstates(1,7);
vx7 = pstates(2,7);
vy7 = pstates(3,7);
vz7 = pstates(4,7);
pg7 = pstates(5,7);
by7 = pstates(6,7);
bz7 = pstates(7,7);
bp7 = bperp(7);

V8 = 1/pstates(1,8);
vx8 = pstates(2,8);
vy8 = pstates(3,8);
vz8 = pstates(4,8);
pg8 = pstates(5,8);
by8 = pstates(6,8);
bz8 = pstates(7,8);
bp8 = bperp(8);


dby_d = 1;
dby_u = 1;

%-------------------------------------------------------------
%* DJ(1,1) = dvx(2)/dbperp(2)
%* DJ(2,1) = dpg(2)/dbperp(2) + bperp(2)*dbperp(2)/dbperp(2)
%-------------------------------------------------------------

%region 1
dV1bp2 = 0;
dvx1bp2 = 0;
dvy1bp2 = 0;
dvz1bp2 = 0;
dpg1bp2 = 0;
dby1bp2 = 0;
dbz1bp2 = 0;
dbp1bp2 = 0;

% compute fast wave speed derivative w.r.t. bperp(2) between regions 1 and 2	  
V_u = V1;
pg_u = pg1;
by_u = by1;
bz_u = bz1;
bperp_u = bp1;
dV_u = dV1bp2;
dpg_u = dpg1bp2;
dbperp_u = dbp1bp2;

bperp_d = bp2;

dbp2bp2 = 1;
dbperp_d = dbp2bp2;

wfs = wmsl;

Lagrangian_shock_speed_derivative;
dwmslbp2 = dwf;

%* fprintf('\n')
%* fprintf('1. RJ: dwflbp2 = %f\n',dwflbp2);

% Region 2

dby2bp2 = cos(psi(1));

dbz2bp2 = sin(psi(1));

if pg2 <= pg1
  % fast rarefaction
  R0 = pg1*V1^gamma;
  R1 = bp1*V1;
  K1 = sqrt(gamma*R0);
  K2 = (R1/K1)^2;
  
  dvx2bp2 = -K1*bp2^(-2/3)*sqrt(1 + K2*(bp2/(V1*bp1))^(1/3))/(V1*bp1)^(1/3);
  
  dpg2bp2 = gamma*pg1*bp2^(gamma - 1)/bp1^gamma;
  
else  
  % fast shock
  
  dvy2bp2 = dvy1bp2 + (bx/(wmsl*wmsl))*(by2 - by1)*dwmslbp2 ...
	  - (bx/wmsl)*(dby2bp2 - dby1bp2);
	  
  dvz2bp2 = dvz1bp2 + (bx/(wmsl*wmsl))*(bz2 - bz1)*dwmslbp2 ...
	  - (bx/wmsl)*(dbz2bp2 - dbz1bp2);
	  
  dV2bp2 = ((bx*(vy2-vy1)/(wmsl*by2*by2)) - (by1*V1/(by2*by2)))*dby2bp2  ... 
	  + (by1/by2)*dV1bp2 + (V1/by2)*dby1bp2 ... 
          + (bx/(by2*wmsl*wmsl))*(vy2 - vy1)*dwmslbp2 ...
	  - (bx/(by2*wmsl))*(dvy2bp2 - dvy1bp2);
  
%*   dV2bp2 = V1*(dby1bp2 + dbz1bp2)/(by2 + bz2) ...
%* 	  + (by1 + bz1)*dV1bp2/(by2 + bz2) ...
%* 	  - V1*(by1 + bz1)*(dby2bp2 + dbz2bp2)/(by2 + bz2)^2 ...
%* 	  - bx*(dvy2bp2 + dvz2bp2 - dvy1bp2 - dvz1bp2)/(wmsl*(by2 + bz2)) ...
%* 	  + bx*(vy2 + vz2 - vy1 - vz1)*(dwmslbp2/wmsl ... 
%* 			+ (dby2bp2 + dbz2bp2)/(by2 + bz2))/(wmsl*(by2 + bz2));
	 
  dvx2bp2 = dvx1bp2 - (V2 - V1)*dwmslbp2 - wmsl*(dV2bp2 - dV1bp2);

  dpg2bp2 = dpg1bp2 - (bp2*dbp2bp2 - bp1*dbp1bp2) ...
	  + (vx2 - vx1)*dwmslbp2 ...
	  + wmsl*(dvx2bp2 - dvx1bp2);
	  	  
end
	  
DJ(1,1) = dvx2bp2;
DJ(2,1) = dpg2bp2 + bp2*dbp2bp2;

%--------------------------------------------------------------
%* DJ(1,2) = -dvx(7)/dbperp(7)
%* DJ(2,2) = -dpg(7)/dbperp(7) - bperp(7)*dbperp(7)/dbperp(7) 
%--------------------------------------------------------------

%region 8
dV8bp7 = 0;
dvx8bp7 = 0;
dvy8bp7 = 0;
dvz8bp7 = 0;
dpg8bp7 = 0;
dby8bp7 = 0;
dbz8bp7 = 0;
dbp8bp7 = 0;

% compute fast wave speed derivative w.r.t. bperp(2) between regions 7 and 8	  
V_u = V8;
pg_u = pg8;
by_u = by8;
bz_u = bz8;
bperp_u = bp8;
dV_u = dV8bp7;
dpg_u = dpg8bp7;
dbperp_u = dbp8bp7;

bperp_d = bp7;

dbp7bp7 = 1;
dbperp_d = dbp7bp7;

wfs = wmsr;

Lagrangian_shock_speed_derivative;
dwmsrbp7 = dwf;

% Region 7
dby7bp7 = cos(psi(end));

dbz7bp7 = sin(psi(end));

if pg7 <= pg8
  % fast rarefaction
  R0 = pg8*V8^gamma;
  R1 = bp8*V8;
  K1 = sqrt(gamma*R0);
  K2 = (R1/K1)^2;
  
  dvx7bp7 = K1*bp7^(-2/3)*sqrt(1 + K2*(bp7/(V8*bp8))^(1/3))/(V8*bp8)^(1/3);
  
  dpg7bp7 = gamma*pg8*bp7^(gamma - 1)/bp8^gamma;
  
else  

  dvy7bp7 = dvy8bp7 + (bx/(wmsr*wmsr))*(by7 - by8)*dwmsrbp7 ...
	  - (bx/wmsr)*(dby7bp7 - dby8bp7);
	  
  dvz7bp7 = dvz8bp7 + (bx/(wmsr*wmsr))*(bz7 - bz8)*dwmsrbp7 ...
	  - (bx/wmsr)*(dbz7bp7 - dbz8bp7);
	  
  dV7bp7 = ((bx*(vy7 - vy8)/(wmsr*by7*by7)) - (by8*V8/(by7*by7)))*dby7bp7  ... 
	  + (by8/by7)*dV8bp7 + (V8/by7)*dby8bp7 ... 
	 + (bx/(by7*wmsr*wmsr))*(vy7 - vy8)*dwmsrbp7 ...
	 - (bx/(by7*wmsr))*(dvy7bp7 - dvy8bp7);	 
	  
%*   dV7bp7 = ((bx*(vz7 - vz8)/(wmsr*bz7*bz7)) - (bz8*V8/(bz7*bz7)))*dbz7bp7 ... 
%* 	  + (bz8/bz7)*dV8bp7 + (V8/bz7)*dbz8bp7 ... 
%* 	 + (bx/(bz7*wmsr*wmsr))*(vz7 - vz8)*dwmsrbp7 ...
%* 	 - (bx/(bz7*wmsr))*(dvz7bp7 - dvz8bp7);	 	  
%* 	  	  
%*   dV7bp7 = V8*(dby8bp7 + dbz8bp7)/(by7 + bz7) ...
%* 	  + (by8 + bz8)*dV8bp7/(by7 + bz7) ...
%* 	  - V8*(by8 + bz8)*(dby7bp7 + dbz7bp7)/(by7 + bz7)^2 ...
%* 	  - bx*(dvy7bp7 + dvz7bp7 - dvy8bp7 - dvz8bp7)/(wmsr*(by7 + bz7)) ...
%* 	  + bx*(vy7 + vz7 - vy8 - vz8)*(dwmsrbp7/wmsr ... 
%* 			+ (dby7bp7 + dbz7bp7)/(by7 + bz7))/(wmsr*(by7 + bz7));	  	  
	 
  dvx7bp7 = dvx8bp7 - (V7 - V8)*dwmsrbp7 - wmsr*(dV7bp7 - dV8bp7);

  dpg7bp7 = dpg8bp7 - (bp7*dbp7bp7 - bp8*dbp8bp7) ...
	  + (vx7 - vx8)*dwmsrbp7 ...
	  + wmsr*(dvx7bp7 - dvx8bp7);

end

DJ(1,2) = -dvx7bp7;
DJ(2,2) = -dpg7bp7 - bp7*dbp7bp7;












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if 1 == 0
  
%region 4
dV4bp3 = 0;
dvx4bp3 = 0;
dvy4bp3 = 0;
dvz4bp3 = 0;
dpg4bp3 = 0;
dby4bp3 = 0;
dbz4bp3 = 0;
dbp4bp3 = 0;

% compute fast wave speed derivative w.r.t. bperp(3) between regions 3 and 4	  
V_u = V4;
pg_u = pg4;
by_u = by4;
bz_u = bz4;
bperp_u = bp4;
dV_u = dV4bp3;
dpg_u = dpg4bp3;
dbperp_u = dbp4bp3;

bperp_d = bp3;

dbp3bp3 = 1;
dbperp_d = dbp3bp3;

wfs = wmsr;

Lagrangian_shock_speed_derivative;
dwmsrbp3 = dwf;

% Region 3
dby3bp3 = cos(psi(end));

dbz3bp3 = sin(psi(end));

if pg3 <= pg4
  % fast rarefaction
  R0 = pg4*V4^gamma;
  R1 = bp4*V4;
  K1 = sqrt(gamma*R0);
  K2 = (R1/K1)^2;
  
  dvx3bp3 = K1*bp3^(-2/3)*sqrt(1 + K2*(bp3/(V4*bp4))^(1/3))/(V4*bp4)^(1/3);
  
  dpg3bp3 = gamma*pg4*bp3^(gamma - 1)/bp4^gamma;
  
else  
  % fast shock
  dvy3bp3 = dvy4bp3 + (bx/(wmsr*wmsr))*(by3 - by4)*dwmsrbp3 ...
	  - (bx/wmsr)*(dby3bp3 - dby4bp3);
	  
  dvz3bp3 = dvz4bp3 + (bx/(wmsr*wmsr))*(bz3 - bz4)*dwmsrbp3 ...
	  - (bx/wmsr)*(dbz3bp3 - dbz4bp3);
	  
  dV3bp3 = ((bx*(vy3 - vy4)/(wmsr*by3*by3)) - (by4*V4/(by3*by3)))*dby3bp3  ... 
	  + (by4/by3)*dV4bp3 + (V4/by3)*dby4bp3 ... 
	 + (bx/(by3*wmsr*wmsr))*(vy3 - vy4)*dwmsrbp3 ...
	 - (bx/(by3*wmsr))*(dvy3bp3 - dvy4bp3);	 
  
%*   dV3bp3 = V4*(dby4bp3 + dbz4bp3)/(by3 + bz3) ...
%* 	  + (by4 + bz4)*dV4bp3/(by3 + bz3) ...
%* 	  - V4*(by4 + bz4)*(dby3bp3 + dbz3bp3)/(by3 + bz3)^2 ...
%* 	  - bx*(dvy3bp3 + dvz3bp3 - dvy4bp3 - dvz4bp3)/(wmsr*(by3 + bz3)) ...
%* 	  + bx*(vy3 + vz3 - vy4 - vz4)*(dwmsrbp3/wmsr ... 
%* 			+ (dby3bp3 + dbz3bp3)/(by3 + bz3))/(wmsr*(by3 + bz3));	    
	 
  dvx3bp3 = dvx4bp3 - (V3 - V4)*dwmsrbp3 - wmsr*(dV3bp3 - dV4bp3);

  dpg3bp3 = dpg4bp3 - (bp3*dbp3bp3 - bp4*dbp4bp3) ...
	  + (vx3 - vx4)*dwmsrbp3 ...
	  + wmsr*(dvx3bp3 - dvx4bp3);

end

DJ(1,2) = -dvx3bp3;
DJ(2,2) = -dpg3bp3 - bp3*dbp3bp3;

end