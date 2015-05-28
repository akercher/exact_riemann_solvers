%******************************************************************************* 
%* Program: fs_jacob.m
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
%* Calculates the jacobian matrix for the system given by equations (11a)
%* - (11d) in [1].  Following [2] the y-component of the magnetic field
%* is replaced with the perpendicular component.
%* 
%* Shock speeds: calculated in rstates_mhd.m
%*     wfl = fast wave moving in negative direction (left).
%*    wrdl = rotational discontinuity (Alven wave) moving in negative
%*           direction (left).
%*     wsl = slow wave moving in negative direction (left).
%*     wsr = slow wave moving in positive direction (right).
%*    wrdr = rotational discontinuity (Alven wave) moving in positive
%*           direction (right).
%*     wfr = fast wave moving in positive direction (right).
%*          
%******************************************************************************* 

% initialize jacobian
DJ = zeros(4,4);

% define state variables in regions 1,2,3,4,5,6,7,8.
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

V3 = 1/pstates(1,3);
vx3 = pstates(2,3);
vy3 = pstates(3,3);
vz3 = pstates(4,3);
pg3 = pstates(5,3);
by3 = pstates(6,3);
bz3 = pstates(7,3);
bp3 = bperp(3);

V4 = 1/pstates(1,4);
vx4 = pstates(2,4);
vy4 = pstates(3,4);
vz4 = pstates(4,4);
pg4 = pstates(5,4);
by4 = pstates(6,4);
bz4 = pstates(7,4);
bp4 = bperp(4);

V5 = 1/pstates(1,5);
vx5 = pstates(2,5);
vy5 = pstates(3,5);
vz5 = pstates(4,5);
pg5 = pstates(5,5);
by5 = pstates(6,5);
bz5 = pstates(7,5);
bp5 = bperp(5);

V6 = 1/pstates(1,6);
vx6 = pstates(2,6);
vy6 = pstates(3,6);
vz6 = pstates(4,6);
pg6 = pstates(5,6);
by6 = pstates(6,6);
bz6 = pstates(7,6);
bp6 = bperp(6);

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


%* dwslbp2 = 1;
%* dwalbp2 = 1;
%* dwflbp2 = 1;
%* 
%* dwslbp4 = 1;
%* 
%* dwsrbp4 = 1;
%* 
%* dwsrbp7 = 1;
%* dwarbp7 = 1;
%* dwfrbp7 = 1;
%* 
%* dwslpsi = 1;
%* dwflpsi = 1;
%* 
%* dwsrpsi = 1;
%* dwfrpsi = 1;

dby_d = 1;
dby_u = 1;

%-------------------------------------------------------------
%* DJ(1,1) = dvx(4)/dbperp(2)
%* DJ(1,2) = dvy(4)/dbperp(2)
%* DJ(1,3) = dvz(4)/dbperp(2)
%* DJ(1,4) = dpg(4)/dbperp(2) + bperp(4)*dbperp(4)/dbperp(2)
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

wfs = wfl;
Lagrangian_shock_speed_derivative;
dwflbp2 = dwf;

% Region 2

dby2bp2 = cos(psi(1));

dbz2bp2 = sin(psi(1));

dvy2bp2 = dvy1bp2 + (bx/(wfl*wfl))*(by2 - by1)*dwflbp2 ...
	  - (bx/wfl)*(dby2bp2 - dby1bp2);
	  
dvz2bp2 = dvz1bp2 + (bx/(wfl*wfl))*(bz2 - bz1)*dwflbp2 ...
	  - (bx/wfl)*(dbz2bp2 - dbz1bp2);
	  
%* dV2bp2 = ((bx*(vy2 - vy1)/(wfl*by2*by2)) - (by1*V1/(by2*by2)))*dby2bp2  ... 
%* 	  + (by1/by2)*dV1bp2 + (V1/by2)*dby1bp2 ... 
%* 	 + (bx/(by2*wfl*wfl))*(vy2 - vy1)*dwflbp2 ...
%* 	 - (bx/(by2*wfl))*(dvy2bp2 - dvy1bp2);	 
%* 	  
%* dV2bp2 = ((bx*(vz2 - vz1)/(wfl*bz2*bz2)) - (bz1*V1/(bz2*bz2)))*dbz2bp2  ... 
%* 	  + (bz1/bz2)*dV1bp2 + (V1/bz2)*dbz1bp2 ... 
%* 	 + (bx/(bz2*wfl*wfl))*(vz2 - vz1)*dwflbp2 ...
%* 	 - (bx/(bz2*wfl))*(dvz2bp2 - dvz1bp2);	 
	 
dV2bp2 = V1*(dby1bp2 + dbz1bp2)/(by2 + bz2) ...
	  + (by1 + bz1)*dV1bp2/(by2 + bz2) ...
	  - V1*(by1 + bz1)*(dby2bp2 + dbz2bp2)/(by2 + bz2)^2 ...
	  - bx*(dvy2bp2 + dvz2bp2 - dvy1bp2 - dvz1bp2)/(wfl*(by2 + bz2)) ...
	  + bx*(vy2 + vz2 - vy1 - vz1)*(dwflbp2/wfl ... 
			+ (dby2bp2 + dbz2bp2)/(by2 + bz2))/(wfl*(by2 + bz2));
	  
dvx2bp2 = dvx1bp2 - (V2 - V1)*dwflbp2 - wfl*(dV2bp2 - dV1bp2);

dpg2bp2 = dpg1bp2 - (bp2*dbp2bp2 - bp1*dbp1bp2) ...
	  + (vx2 - vx1)*dwflbp2 ...
	  + wfl*(dvx2bp2 - dvx1bp2);

% compute Alfven wave speed derivative w.r.t. bperp(2) between regions 2 and 3	  	  
Cal = sqrt(bx*bx/V2);
Cal = -Cal;
dwalbp2 = -Cal/(2*V2)*dV2bp2;

% Region 3
dbp3bp2 = 1;

dby3bp2 = cos(psi(3));

dbz3bp2 = sin(psi(3));

dV3bp2 = dV2bp2;

dvx3bp2 = dvx2bp2;

dpg3bp2 = dpg2bp2;

dvy3bp2 = dvy2bp2 + (bx/(wal*wal))*(by3 - by2)*dwalbp2 ...
	  - (bx/wal)*(dby3bp2 - dby2bp2);
	  
dvz3bp2 = dvz2bp2 + (bx/(wal*wal))*(bz3 - bz2)*dwalbp2 ...
	  - (bx/wal)*(dbz3bp2 - dbz2bp2);

% compute slow wave speed derivative w.r.t. bperp(2) between regions 3 and 4	  
V_u = V3;
pg_u = pg3;
by_u = by3;
bz_u = bz3;
bperp_u = bp3;
dV_u = dV3bp2;
dpg_u = dpg3bp2;
dbperp_u = dbp3bp2;

bperp_d = bp4;

dbp4bp2 = 0;
dbperp_d = dbp4bp2;

wfs = wsl;
	  
Lagrangian_shock_speed_derivative;
dwslbp2 = dws;

dC023 = dC02;
dCa23 = dCa2;
dCperp23 = dCperp2;
dCf23 = dCf2;
dCs23 = dCs2;
%* fprintf('3. RJ: dwslbp2 = %f\n',dwslbp2);	  

% Region 4
dby4bp2 = 0;
dbz4bp2 = 0;

%* if pg4 <= pg3
%*   % fast rarefaction
%*   R0 = pg3*V3^gamma;
%*   
%*   % calculate derivatives of speeds squared
%*   [dC024,dCa24,dCperp24,dCs24,dCf24] = wave_speed_derivative(gamma,V4,pg4,bp4,bx,...
%* 						      dV4bp2,dpg4bp2,dbp4bp2,...
%* 						      C024,Ca24,Cperp24);
%*   
%*   
%*   dvx2bp2 = -K1*bp2^(-2/3)*sqrt(1 + K2*(bp2/(V1*bp1))^(1/3))/(V1*bp1)^(1/3);
%*   
%*   dpg2bp2 = gamma*pg1*bp2^(gamma - 1)/bp1^gamma;
%*   
%* else  
  % fast shock
  dvy4bp2 = dvy3bp2 + (bx/(wsl*wsl))*(by4 - by3)*dwslbp2 ...
	  - (bx/wsl)*(dby4bp2 - dby3bp2);
	  
  dvz4bp2 = dvz3bp2 + (bx/(wsl*wsl))*(bz4 - bz3)*dwslbp2 ...
	  - (bx/wsl)*(dbz4bp2 - dbz3bp2);
	  
%*   dV4bp2 = ((bx*(vy4-vy3)/(wsl*by4*by4)) - (by3*V3/(by4*by4)))*dby4bp2 ... 
%* 	  + (by3/by4)*dV3bp2 + (V3/by4)*dby3bp2 ... 
%* 	 + (bx/(by4*wsl*wsl))*(vy4 - vy3)*dwslbp2 ...
%* 	 - (bx/(by4*wsl))*(dvy4bp2 - dvy3bp2);	 
%* 
%*   dV4bp2 = ((bx*(vz4-vz3)/(wsl*bz4*bz4)) - (by3*V3/(bz4*bz4)))*dbz4bp2 ... 
%* 	  + (bz3/bz4)*dV3bp2 + (V3/bz4)*dbz3bp2 ... 
%* 	 + (bx/(bz4*wsl*wsl))*(vz4 - vz3)*dwslbp2 ...
%* 	 - (bx/(bz4*wsl))*(dvz4bp2 - dvz3bp2);	 

  dV4bp2 = V3*(dby3bp2 + dbz3bp2)/(by4 + bz4) ...
	  + (by3 + bz3)*dV3bp2/(by4 + bz4) ...
	  - V3*(by3 + bz3)*(dby4bp2 + dbz4bp2)/(by4 + bz4)^2 ...
	  - bx*(dvy4bp2 + dvz4bp2 - dvy3bp2 - dvz3bp2)/(wsl*(by4 + bz4)) ...
	  + bx*(vy4 + vz4 - vy3 - vz3)*(dwslbp2/wsl ... 
			+ (dby4bp2 + dbz4bp2)/(by4 + bz4))/(wsl*(by4 + bz4));
	 
  dvx4bp2 = dvx3bp2 - (V4 - V3)*dwslbp2 - wsl*(dV4bp2 - dV3bp2);

  dpg4bp2 = dpg3bp2 - (bp4*dbp4bp2 - bp3*dbp3bp2) ...
	  + (vx4 - vx3)*dwslbp2 ...
	  + wsl*(dvx4bp2 - dvx3bp2);
	  
%* end

DJ(1,1) = dvx4bp2;
DJ(2,1) = dvy4bp2;
DJ(3,1) = dvz4bp2;
DJ(4,1) = dpg4bp2 + bp4*dbp4bp2;

%-------------------------------------------------------------
%* DJ(1,2) = dvx(4)/dbperp(4) - dvx(5)/dbperp(4)
%* DJ(2,2) = dvy(4)/dbperp(4) - dvy(5)/dbperp(4)
%* DJ(3,2) = dvz(4)/dbperp(4) - dvz(5)/dbperp(4)
%* DJ(4,2) = dpg(4)/dbperp(4) - dpg(5)/dbperp(4)
%*           + bperp(4)*dbperp(4)/dbperp(4) 
%*           - bperp(5)*dbperp(5)/dbperp(4) 
%-------------------------------------------------------------

% Region 3
dV3bp4 = 0;
dvx3bp4 = 0;
dvy3bp4 = 0;
dvz3bp4 = 0;
dby3bp4 = 0;
dbz3bp4 = 0;
dpg3bp4 = 0;
dbp3bp4 = 0;

% compute slow wave speed derivative w.r.t. bperp(4) between regions 3 and 4
V_u = V3;
pg_u = pg3;
by_u = by3;
bz_u = bz3;
bperp_u = bp3;
dV_u = dV3bp4;
dpg_u = dpg3bp4;
dbperp_u = dbp3bp4;

bperp_d = bp4;
dbp4bp4 = 1;
dbperp_d = dbp4bp4;

wfs = wsl;
	  
Lagrangian_shock_speed_derivative;
dwslbp4 = dws;

% Region 4
dby4bp4 = cos(psi(3));

dbz4bp4 = sin(psi(3));

dvy4bp4 = dvy3bp4 + (bx/(wsl*wsl))*(by4 - by3)*dwslbp4 ...
	  - (bx/wsl)*(dby4bp4 - dby3bp4);
	  
dvz4bp4 = dvz3bp4 + (bx/(wsl*wsl))*(bz4 - bz3)*dwslbp4 ...
	  - (bx/wsl)*(dbz4bp4 - dbz3bp4);
	  
%* dV4bp4 = ((bx*(vy4 - vy3)/(wsl*by4*by4)) - (by3*V3/(by4*by4)))*dby4bp4 ... 
%* 	  + (by3/by4)*dV3bp4 + (V3/by4)*dby3bp4 ... 
%* 	 + (bx/(by4*wsl*wsl))*(vy4 - vy3)*dwslbp4 ...
%* 	 - (bx/(by4*wsl))*(dvy4bp4 - dvy3bp4);	 
%* 	  
%* dV4bp4 = ((bx*(vz4 - vz3)/(wsl*bz4*bz4)) - (bz3*V3/(bz4*bz4)))*dbz4bp4 ... 
%* 	  + (bz3/bz4)*dV3bp4 + (V3/bz4)*dbz3bp4 ... 
%* 	 + (bx/(bz4*wsl*wsl))*(vz4 - vz3)*dwslbp4 ...
%* 	 - (bx/(bz4*wsl))*(dvz4bp4 - dvz3bp4);
	  
dV4bp4 = V3*(dby3bp4 + dbz3bp4)/(by4 + bz4) ...
	  + (by3 + bz3)*dV3bp4/(by4 + bz4) ...
	  - V3*(by3 + bz3)*(dby4bp4 + dbz4bp4)/(by4 + bz4)^2 ...
	  - bx*(dvy4bp4 + dvz4bp4 - dvy3bp4 - dvz3bp4)/(wsl*(by4 + bz4)) ...
	  + bx*(vy4 + vz4 - vy3 - vz3)*(dwslbp4/wsl ... 
			+ (dby4bp4 + dbz4bp4)/(by4 + bz4))/(wsl*(by4 + bz4));	  
	  
dvx4bp4 = dvx3bp4 - (V4 - V3)*dwslbp4 - wsl*(dV4bp4 - dV3bp4);

dpg4bp4 = dpg3bp4 - (bp4*dbp4bp4 - bp3*dbp3bp4) ...
	  + (vx4 - vx3)*dwslbp4 ...
	  + wsl*(dvx4bp4 - dvx3bp4);

% Region 6
dV6bp4 = 0;
dvx6bp4 = 0;
dvy6bp4 = 0;
dvz6bp4 = 0;	  
dby6bp4 = 0;
dbz6bp4 = 0;
dpg6bp4 = 0;
dbp6bp4 = 0;

% compute slow wave speed derivative w.r.t. bperp(4) between regions 5 and 6
V_u = V6;
pg_u = pg6;
by_u = by6;
bz_u = bz6;
bperp_u = bp6;
dV_u = dV6bp4;
dpg_u = dpg6bp4;
dbperp_u = dbp6bp4;

bperp_d = bp5;
dbp5bp4 = 1;
dbperp_d = dbp5bp4;

wfs = wsr;

Lagrangian_shock_speed_derivative;
dwsrbp4 = dws;

% Region 5
dby5bp4 = cos(psi(3));

dbz5bp4 = sin(psi(3));

dvy5bp4 = dvy6bp4 + (bx/(wsr*wsr))*(by5 - by6)*dwsrbp4 ...
	  - (bx/wsr)*(dby5bp4 - dby6bp4);
	  
dvz5bp4 = dvz6bp4 + (bx/(wsr*wsr))*(bz5 - bz6)*dwsrbp4 ...
	  - (bx/wsr)*(dbz5bp4 - dbz6bp4);
	  
%* dV5bp4 = ((bx*(vy5 - vy6)/(wsr*by5*by5)) - (by6*V6/(by5*by5)))*dby5bp4 ... 
%* 	  + (by6/by5)*dV6bp4 + (V6/by5)*dby6bp4 ... 
%* 	 + (bx/(by5*wsr*wsr))*(vy5 - vy6)*dwsrbp4 ...
%* 	 - (bx/(by5*wsr))*(dvy5bp4 - dvy6bp4);	 
	  
dV5bp4 = V6*(dby6bp4 + dbz6bp4)/(by5 + bz5) ...
	  + (by6 + bz6)*dV6bp4/(by5 + bz5) ...
	  - V6*(by6 + bz6)*(dby5bp4 + dbz5bp4)/(by5 + bz5)^2 ...
	  - bx*(dvy5bp4 + dvz5bp4 - dvy6bp4 - dvz6bp4)/(wsr*(by5 + bz5)) ...
	  + bx*(vy5 + vz5 - vy6 - vz6)*(dwsrbp4/wsr ... 
			+ (dby5bp4 + dbz5bp4)/(by5 + bz5))/(wsr*(by5 + bz5));	  
	 
dvx5bp4 = dvx6bp4 - (V5 - V6)*dwsrbp4 - wsr*(dV5bp4 - dV6bp4);

dpg5bp4 = dpg6bp4 - (bp5*dbp5bp4 - bp6*dbp6bp4) ...
	  + (vx5 - vx6)*dwsrbp4 ...
	  + wsr*(dvx5bp4 - dvx6bp4);
	  
DJ(1,2) = dvx4bp4 - dvx5bp4;
DJ(2,2) = dvy4bp4 - dvy5bp4;
DJ(3,2) = dvz4bp4 - dvz5bp4;
DJ(4,2) = dpg4bp4 - dpg5bp4 + bp4*dbp4bp4 - bp5*dbp5bp4;

%--------------------------------------------------------------
%* DJ(1,3) = -dvx(5)/dbperp(7)
%* DJ(2,3) = -dvy(5)/dbperp(7)
%* DJ(3,3) = -dvz(5)/dbperp(7)
%* DJ(4,3) = -dpg(5)/dbperp(7) - bperp(5)*dbperp(5)/dbperp(7) 
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

% compute fast wave speed derivative w.r.t. bperp(7) between regions 7 and 8	  
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

wfs = wfr;

Lagrangian_shock_speed_derivative;
dwfrbp7 = dwf;

% Region 7
dby7bp7 = cos(psi(end));

dbz7bp7 = sin(psi(end));

dvy7bp7 = dvy8bp7 + (bx/(wfr*wfr))*(by7 - by8)*dwfrbp7 ...
	  - (bx/wfr)*(dby7bp7 - dby8bp7);
	  
dvz7bp7 = dvz8bp7 + (bx/(wfr*wfr))*(bz7 - bz8)*dwfrbp7 ...
	  - (bx/wfr)*(dbz7bp7 - dbz8bp7);
	  
%* dV7bp7 = ((bx*(vy7 - vy8)/(wfr*by7*by7)) - (by8*V8/(by7*by7)))*dby7bp7  ... 
%* 	  + (by8/by7)*dV8bp7 + (V8/by7)*dby8bp7 ... 
%* 	 + (bx/(by7*wfr*wfr))*(vy7 - vy8)*dwfrbp7 ...
%* 	 - (bx/(by7*wfr))*(dvy7bp7 - dvy8bp7);	 
%* 	  
%* dV7bp7 = ((bx*(vz7 - vz8)/(wfr*bz7*bz7)) - (bz8*V8/(bz7*bz7)))*dbz7bp7 ... 
%* 	  + (bz8/bz7)*dV8bp7 + (V8/bz7)*dbz8bp7 ... 
%* 	 + (bx/(bz7*wfr*wfr))*(vz7 - vz8)*dwfrbp7 ...
%* 	 - (bx/(bz7*wfr))*(dvz7bp7 - dvz8bp7);	 	  
	  	  
dV7bp7 = V8*(dby8bp7 + dbz8bp7)/(by7 + bz7) ...
	  + (by8 + bz8)*dV8bp7/(by7 + bz7) ...
	  - V8*(by8 + bz8)*(dby7bp7 + dbz7bp7)/(by7 + bz7)^2 ...
	  - bx*(dvy7bp7 + dvz7bp7 - dvy8bp7 - dvz8bp7)/(wfr*(by7 + bz7)) ...
	  + bx*(vy7 + vz7 - vy8 - vz8)*(dwfrbp7/wfr ... 
			+ (dby7bp7 + dbz7bp7)/(by7 + bz7))/(wfr*(by7 + bz7));	  	  
	 
dvx7bp7 = dvx8bp7 - (V7 - V8)*dwfrbp7 - wfr*(dV7bp7 - dV8bp7);

dpg7bp7 = dpg8bp7 - (bp7*dbp7bp7 - bp8*dbp8bp7) ...
	  + (vx7 - vx8)*dwfrbp7 ...
	  + wfr*(dvx7bp7 - dvx8bp7);


% compute Alfven wave speed derivative w.r.t. bperp(7) between regions 6 and 7
Car = sqrt(bx*bx/V7);
dwarbp7 = -Car/(2*V7)*dV7bp7;
	  
% Region 6
dbp6bp7 = 1;

dby6bp7 = cos(psi(3));

dbz6bp7 = sin(psi(3));

dV6bp7 = dV7bp7;

dvx6bp7 = dvx7bp7;

dpg6bp7 = dpg7bp7;

dvy6bp7 = dvy7bp7 + (bx/(war*war))*(by6 - by7)*dwarbp7 ...
	  - (bx/war)*(dby6bp7 - dby7bp7);
	  
dvz6bp7 = dvz7bp7 + (bx/(war*war))*(bz6 - bz7)*dwarbp7 ...
	  - (bx/war)*(dbz6bp7 - dbz7bp7);

% compute slow wave speed derivative w.r.t. bperp(7) between regions 5 and 6	  
V_u = V6;
pg_u = pg6;
by_u = by6;
bz_u = bz6;
bperp_u = bp6;
dV_u = dV6bp7;
dpg_u = dpg6bp7;
dbperp_u = dbp6bp7;

bperp_d = bp5;
dbp5bp7 = 0;
dbperp_d = dbp5bp7;

wfs = wsr;
	  
Lagrangian_shock_speed_derivative;
dwsrbp7 = dws;
	  
% Region 5
dby5bp7 = 0;
dbz5bp7 = 0;

dvy5bp7 = dvy6bp7 + (bx/(wsr*wsr))*(by5 - by6)*dwsrbp7 ...
	  - (bx/wsr)*(dby5bp7 - dby6bp7);
	  
dvz5bp7 = dvz6bp7 + (bx/(wsr*wsr))*(bz5 - bz6)*dwsrbp7 ...
	  - (bx/wsr)*(dbz5bp7 - dbz6bp7);

%* dV5bp7 = ((bx*(vy5 - vy6)/(wsr*by5*by5)) - (by6*V6/(by5*by5)))*dby5bp7 ... 
%* 	  + (by6/by5)*dV6bp7 + (V6/by5)*dby6bp7 ... 
%* 	 + (bx/(by5*wsr*wsr))*(vy5 - vy6)*dwsrbp7 ...
%* 	 - (bx/(by5*wsr))*(dvy5bp7 - dvy6bp7);	 	  
%* 	  
%* dV5bp7 = ((bx*(vz5 - vz6)/(wsr*bz5*bz5)) - (bz6*V6/(bz5*bz5)))*dbz5bp7 ... 
%* 	  + (bz6/bz5)*dV6bp7 + (V6/bz5)*dbz6bp7 ... 
%* 	 + (bx/(bz5*wsr*wsr))*(vz5 - vz6)*dwsrbp7 ...
%* 	 - (bx/(bz5*wsr))*(dvz5bp7 - dvz6bp7);	 	  
	  	  
dV5bp7 = V6*(dby6bp7 + dbz6bp7)/(by5 + bz5) ...
	  + (by6 + bz6)*dV6bp7/(by5 + bz5) ...
	  - V6*(by6 + bz6)*(dby5bp7 + dbz5bp7)/(by5 + bz5)^2 ...
	  - bx*(dvy5bp7 + dvz5bp7 - dvy6bp7 - dvz6bp7)/(wsr*(by5 + bz5)) ...
	  + bx*(vy5 + vz5 - vy6 - vz6)*(dwsrbp7/wsr ... 
			+ (dby5bp7 + dbz5bp7)/(by5 + bz5))/(wsr*(by5 + bz5));	  
	  	 
dvx5bp7 = dvx6bp7 - (V5 - V6)*dwsrbp7 - wsr*(dV5bp7 - dV6bp7);

dpg5bp7 = dpg6bp7 - (bp5*dbp5bp7 - bp6*dbp6bp7) ...
	  + (vx5 - vx6)*dwsrbp7 ...
	  + wsr*(dvx5bp7 - dvx6bp7);

DJ(1,3) = -dvx5bp7;
DJ(2,3) = -dvy5bp7;
DJ(3,3) = -dvz5bp7;
DJ(4,3) = -dpg5bp7 - bp5*dbp5bp7;

%-------------------------------------------------------------
%* DJ(1,4) = dvx(4)/dpsi - dvx(5)/dpsi
%* DJ(2,4) = dvy(4)/dpsi - dvy(5)/dpsi
%* DJ(3,4) = dvz(4)/dpsi - dvz(5)/dpsi
%* DJ(4,4) = dpg(4)/dpsi - dpg(5)/dpsi
%*           + bperp(4)*dbperp(4)/dpsi
%*           - bperp(5)*dbperp(5)/dpsi 
%-------------------------------------------------------------

%region 2
dV2psi = 0;
dvx2psi = 0;
dvy2psi = 0;
dvz2psi = 0;
dpg2psi = 0;
dby2psi = 0;
dbz2psi = 0;
dbp2psi = 0;

% compute Alfven wave speed derivative between regions 2 and 3
Cal = sqrt(bx*bx/V2);
Cal = -Cal;
dwalpsi = -Cal/(2*V2)*dV2psi;

% Region 3
dV3psi = 0;
dvx3psi = 0;
dpg3psi = 0;
dbp3psi = 0;

dby3psi = -bp2*sin(psi(3));

dbz3psi = bp2*cos(psi(3));


dvy3psi = dvy2psi + (bx/(wal*wal))*(by3 - by2)*dwalpsi ...
	  - (bx/wal)*(dby3psi - dby2psi);
	  
dvz3psi = dvz2psi + (bx/(wal*wal))*(bz3 - bz2)*dwalpsi ...
	  - (bx/wal)*(dbz3psi - dbz2psi);

% compute slow wave speed derivative w.r.t. psi between regions 3 and 4	  
V_u = V3;
pg_u = pg3;
by_u = by3;
bz_u = bz3;
bperp_u = bp3;
dV_u = dV3psi;
dpg_u = dpg3psi;
dbperp_u = dbp3psi;

bperp_d = bp4;

dbp4psi = 0;
dbperp_d = dbp4psi;

wfs = wsl;
	  
Lagrangian_shock_speed_derivative;
dwslpsi = dws;
	  
% Region 4
dby4psi = -bp4*sin(psi(3));
dbz4psi = bp4*cos(psi(3));

dvy4psi = dvy3psi + (bx/(wsl*wsl))*(by4 - by3)*dwslpsi ...
	  - (bx/wsl)*(dby4psi - dby3psi);
	  
dvz4psi = dvz3psi + (bx/(wsl*wsl))*(bz4 - bz3)*dwslpsi ...
	  - (bx/wsl)*(dbz4psi - dbz3psi);

%* dV4psi = ((bx*(vy4 - vy3)/(wsl*by4*by4)) - (by3*V3/(by4*by4)))*dby4psi ... 
%* 	  + (by3/by4)*dV3psi + (V3/by4)*dby3psi ... 
%* 	 + (bx/(by4*wsl*wsl))*(vy4 - vy3)*dwslpsi ...
%* 	 - (bx/(by4*wsl))*(dvy4psi - dvy3psi);	 
	  
%* dV4psi = ((bx*(vz4 - vz3)/(wsl*bz4*bz4)) - (bz3*V3/(bz4*bz4)))*dbz4psi ... 
%* 	  + (bz3/bz4)*dV3psi + (V3/bz4)*dbz3psi ... 
%* 	 + (bx/(bz4*wsl*wsl))*(vz4 - vz3)*dwslpsi ...
%* 	 - (bx/(bz4*wsl))*(dvz4psi - dvz3psi);	 
	  
dV4psi = V3*(dby3psi + dbz3psi)/(by4 + bz4) ...
	  + (by3 + bz3)*dV3psi/(by4 + bz4) ...
	  - V3*(by3 + bz3)*(dby4psi + dbz4psi)/(by4 + bz4)^2 ...
	  - bx*(dvy4psi + dvz4psi - dvy3psi - dvz3psi)/(wsl*(by4 + bz4)) ...
	  + bx*(vy4 + vz4 - vy3 - vz3)*(dwslpsi/wsl ... 
			+ (dby4psi + dbz4psi)/(by4 + bz4))/(wsl*(by4 + bz4));
	 
dvx4psi = dvx3psi - (V4 - V3)*dwslpsi - wsl*(dV4psi - dV3psi);

dpg4psi = dpg3psi - (bp4*dbp4psi - bp3*dbp3psi) ...
	  + (vx4 - vx3)*dwslpsi ...
	  + wsl*(dvx4psi - dvx3psi);


%region 7
dV7psi = 0;
dvx7psi = 0;
dvy7psi = 0;
dvz7psi = 0;
dpg7psi = 0;
dby7psi = 0;
dbz7psi = 0;
dbp7psi = 0;

% compute Alfven wave speed derivative between regions 6 and 7
Car = sqrt(bx*bx/V7);
dwarpsi = -Car/(2*V7)*dV7psi;

% Region 6
dV6psi = 0;
dvx6psi = 0;
dpg6psi = 0;	  
dbp6psi = 0;

dby6psi = -bp7*sin(psi(3));

dbz6psi = bp7*cos(psi(3));

dvy6psi = dvy7psi + (bx/(war*war))*(by6 - by7)*dwarpsi ...
	  - (bx/war)*(dby6psi - dby7psi);
	  
dvz6psi = dvz7psi + (bx/(war*war))*(bz6 - bz7)*dwarpsi ...
	  - (bx/war)*(dbz6psi - dbz7psi);

% compute slow wave speed derivative w.r.t. psi between regions 5 and 6	  
V_u = V6;
pg_u = pg6;
by_u = by6;
bz_u = bz6;
bperp_u = bp6;
dV_u = dV6psi;
dpg_u = dpg6psi;
dbperp_u = dbp6psi;

bperp_d = bp5;

dbp5psi = 0;
dbperp_d = dbp5psi;

wfs = wsr;

Lagrangian_shock_speed_derivative;
dwsrpsi = dws;
	  
% Region 5
dby5psi = -bp5*sin(psi(3));
dbz5psi = bp5*cos(psi(3));

dvy5psi = dvy6psi + (bx/(wsr*wsr))*(by5 - by6)*dwsrpsi ...
	  - (bx/wsr)*(dby5psi - dby6psi);
	  
dvz5psi = dvz6psi + (bx/(wsr*wsr))*(bz5 - bz6)*dwsrpsi ...
	  - (bx/wsr)*(dbz5psi - dbz6psi);

%* dV5psi = ((bx*(vy5-vy6)/(wsr*by5*by5)) - (by6*V6/(by5*by5)))*dby5psi ... 
%* 	  + (by6/by5)*dV6psi + (V6/by5)*dby5psi ... 
%* 	 + (bx/(by5*wsr*wsr))*(vy5 - vy6)*dwsrpsi ...
%* 	 - (bx/(by5*wsr))*(dvy5psi - dvy6psi);	 
	  
%* dV5psi = ((bx*(vz5-vz6)/(wsr*bz5*bz5)) - (bz6*V6/(bz5*bz5)))*dbz5psi ... 
%* 	  + (bz6/bz5)*dV6psi + (V6/bz5)*dbz5psi ... 
%* 	 + (bx/(bz5*wsr*wsr))*(vz5 - vz6)*dwsrpsi ...
%* 	 - (bx/(bz5*wsr))*(dvz5psi - dvz6psi);	 
	  
dV5psi = V6*(dby6psi + dbz6psi)/(by5 + bz5) ...
	  + (by6 + bz6)*dV6psi/(by5 + bz5) ...
	  - V6*(by6 + bz6)*(dby5psi + dbz5psi)/(by5 + bz5)^2 ...
	  - bx*(dvy5psi + dvz5psi - dvy6psi - dvz6psi)/(wsr*(by5 + bz5)) ...
	  + bx*(vy5 + vz5 - vy6 - vz6)*(dwsrpsi/wsr ... 
			+ (dby5psi + dbz5psi)/(by5 + bz5))/(wsr*(by5 + bz5));	  	  
	 
dvx5psi = dvx6psi - (V5 - V6)*dwsrpsi - wsr*(dV5psi - dV6psi);

dpg5psi = dpg6psi - (bp5*dbp5psi - bp6*dbp6psi) ...
	  + (vx5 - vx6)*dwsrpsi ...
	  + wsr*(dvx5psi - dvx6psi);


DJ(1,4) = dvx4psi - dvx5psi;
DJ(2,4) = dvy4psi - dvy5psi;
DJ(3,4) = dvz4psi - dvz5psi;
DJ(4,4) = dpg4psi - dpg5psi + bp4*dbp4psi - bp5*dbp5psi;

if 1 == 0
invDJ = inv(DJ);
detDJ = det(DJ);

%* fprintf('\n');

b14 = DJ(1,2)*DJ(2,4)*DJ(3,3) + DJ(1,3)*DJ(2,2)*DJ(3,4) ...
      + DJ(1,4)*DJ(2,3)*DJ(3,2) - (DJ(1,2)*DJ(2,3)*DJ(3,4) ...
      + DJ(1,3)*DJ(2,4)*DJ(3,2) + DJ(1,4)*DJ(2,2)*DJ(3,3));

fprintf('\n');

b211 = DJ(2,1)*DJ(3,4)*DJ(4,3)/detDJ
b212 = DJ(2,3)*DJ(3,1)*DJ(4,4)/detDJ
b213 = DJ(2,4)*DJ(3,3)*DJ(4,1)/detDJ 
b214 = DJ(2,1)*DJ(3,3)*DJ(4,4)/detDJ
b215 = DJ(2,3)*DJ(3,4)*DJ(4,1)/detDJ 
b216 = DJ(2,4)*DJ(3,1)*DJ(4,3)/detDJ

b21 = b211 + b212 + b213 - (b214 + b215 + b216);

fprintf('\n');

b221 = DJ(1,1)*DJ(3,3)*DJ(4,4)/detDJ      
b222 = DJ(1,3)*DJ(3,4)*DJ(4,1)/detDJ      
b223 = DJ(1,4)*DJ(3,1)*DJ(4,3)/detDJ      
b224 = DJ(1,1)*DJ(3,4)*DJ(4,3)/detDJ      
b225 = DJ(1,3)*DJ(3,1)*DJ(4,4)/detDJ      
b226 = DJ(1,4)*DJ(3,3)*DJ(4,1)/detDJ      

b22 = b221 + b222 + b223 - (b224 + b225 + b226);

fprintf('\n');

b231 = DJ(1,1)*DJ(2,4)*DJ(4,3)/detDJ
b232 = DJ(1,3)*DJ(2,1)*DJ(4,4)/detDJ      
b233 = DJ(1,4)*DJ(2,3)*DJ(4,1)/detDJ      
b234 = DJ(1,1)*DJ(2,3)*DJ(4,4)/detDJ      
b235 = DJ(1,3)*DJ(2,4)*DJ(4,1)/detDJ      
b236 = DJ(1,4)*DJ(2,1)*DJ(4,3)/detDJ      

b23 = b231 + b232 + b233 - (b234 + b235 + b236);

fprintf('\n');

b24 = DJ(1,1)*DJ(2,2)*DJ(3,4) + DJ(1,3)*DJ(2,4)*DJ(3,1) ...
      + DJ(1,4)*DJ(2,1)*DJ(3,3) - (DJ(1,1)*DJ(2,4)*DJ(3,3) ...
      + DJ(1,3)*DJ(2,1)*DJ(3,4) + DJ(1,4)*DJ(2,3)*DJ(3,1));

b34 = DJ(1,1)*DJ(2,4)*DJ(3,2) + DJ(1,2)*DJ(2,1)*DJ(3,4) ...
      + DJ(1,4)*DJ(2,2)*DJ(3,1) - (DJ(1,1)*DJ(2,2)*DJ(3,4) ...
      + DJ(1,2)*DJ(2,4)*DJ(3,1) + DJ(1,4)*DJ(2,1)*DJ(3,2));

%* b = DJ(,)*DJ(,)*DJ(,)/detDJ      
%* b34 = DJ(,)*DJ(,)*DJ(,) + DJ(,)*DJ(,)*DJ(,) ...
%*       + DJ(,)*DJ(,)*DJ(,) - (DJ(,)*DJ(,)*DJ(,) ...
%*       + DJ(,)*DJ(,)*DJ(,) + DJ(,)*DJ(,)*DJ(,));
%* 
%* b34 = DJ(,)*DJ(,)*DJ(,) + DJ(,)*DJ(,)*DJ(,) ...
%*       + DJ(,)*DJ(,)*DJ(,) - (DJ(,)*DJ(,)*DJ(,) ...
%*       + DJ(,)*DJ(,)*DJ(,) + DJ(,)*DJ(,)*DJ(,));
      
keyboard

end

%* fprintf('\n')
%* fprintf('wfl = %f dwsl/dbp2 = %f\n',wfl,dwflbp2);
%* fprintf('wal = %f dwal/dbp2 = %f dwal/dpsi = %f\n',wal,dwalbp2,dwalpsi);
%* fprintf('wsl = %f dwsl/dbp2 = %f dwsl/dbp4 = %f dwsl/dpsi = %f\n',wsl, ...
%* 	dwslbp2,dwslbp4,dwslpsi)
%* fprintf('wsr = %f dwsr/dbp4 = %f dwsr/dbp7 = %f dwsr/dpsi = %f\n',wsr, ...
%* 	dwsrbp4,dwsrbp7,dwsrpsi)
%* fprintf('war = %f dwar/dbp7 = %f dwar/dpsi = %f\n',war,dwarbp7,dwarpsi);
%* fprintf('wfr = %f dwsr/dbp7 = %f\n',wfr,dwfrbp7);
%* fprintf('\n')









