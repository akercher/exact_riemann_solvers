%******************************************************************************* 
%* Program: fs_jacob_dw.m
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
%* - (11d) in [1].  
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

% define state variables in regions 1-8.
V1 = 1/pstates(1,1);
vx1 = pstates(2,1);
vy1 = pstates(3,1);
vz1 = pstates(4,1);
pg1 = pstates(5,1);
by1 = pstates(6,1);
bz1 = pstates(7,1);
%* bp1 = bperp(1);

V2 = 1/pstates(1,2);
vx2 = pstates(2,2);
vy2 = pstates(3,2);
vz2 = pstates(4,2);
pg2 = pstates(5,2);
by2 = pstates(6,2);
bz2 = pstates(7,2);
%* bp2 = bperp(2);

V3 = 1/pstates(1,3);
vx3 = pstates(2,3);
vy3 = pstates(3,3);
vz3 = pstates(4,3);
pg3 = pstates(5,3);
by3 = pstates(6,3);
bz3 = pstates(7,3);
%* bp3 = bperp(3);

V4 = 1/pstates(1,4);
vx4 = pstates(2,4);
vy4 = pstates(3,4);
vz4 = pstates(4,4);
pg4 = pstates(5,4);
by4 = pstates(6,4);
bz4 = pstates(7,4);
%* bp4 = bperp(4);

V5 = 1/pstates(1,5);
vx5 = pstates(2,5);
vy5 = pstates(3,5);
vz5 = pstates(4,5);
pg5 = pstates(5,5);
by5 = pstates(6,5);
bz5 = pstates(7,5);
%* bp5 = bperp(5);

V6 = 1/pstates(1,6);
vx6 = pstates(2,6);
vy6 = pstates(3,6);
vz6 = pstates(4,6);
pg6 = pstates(5,6);
by6 = pstates(6,6);
bz6 = pstates(7,6);
%* bp6 = bperp(6);

V7 = 1/pstates(1,7);
vx7 = pstates(2,7);
vy7 = pstates(3,7);
vz7 = pstates(4,7);
pg7 = pstates(5,7);
by7 = pstates(6,7);
bz7 = pstates(7,7);
%* bp7 = bperp(7);

V8 = 1/pstates(1,8);
vx8 = pstates(2,8);
vy8 = pstates(3,8);
vz8 = pstates(4,8);
pg8 = pstates(5,8);
by8 = pstates(6,8);
bz8 = pstates(7,8);
%* bp8 = bperp(8);


%* dwslbp2 = 1;
%* dwrdlbp2 = 1;
%* dwflbp2 = 1;
%* 
%* dwslbp4 = 1;
%* 
%* dwsrbp4 = 1;
%* 
%* dwsrbp7 = 1;
%* dwrdrbp7 = 1;
%* dwfrbp7 = 1;
%* 
%* dwslpsi = 1;
%* dwflpsi = 1;
%* 
%* dwsrpsi = 1;
%* dwfrpsi = 1;

%---------------------------------------------------------
%* DJ(1,1) = dvx(4)/dby(2)
%* DJ(1,2) = dvy(4)/dby(2)
%* DJ(1,3) = dvz(4)/dby(2)
%* DJ(1,4) = dpg(4)/dby(2)
%---------------------------------------------------------

%region 1
dV1by2 = 0;
dvx1by2 = 0;
dvy1by2 = 0;
dvz1by2 = 0;
dpg1by2 = 0;
dby1by2 = 0;
dbz1by2 = 0;
dbp1by2 = 0;

% compute fast wave speed derivative w.r.t. by(2) between regions 1 and 2	  
V_d = V1;
pg_d = pg1;
by_d = by1;
bz_d = bz1;
dV_d = dV1by2;
dpg_d = dpg1by2;
dby_d= dby1by2;

dby2by2 = 1;
dby_u= dby2by2;

bperp_d = bperp(1);
bperp_u = bperp(2);
dbperp_d = 0;
dbperp_u = 1;

wfs = wfl;

wave_speed_derivative;
dwflby2 = dwf_dw;

fprintf('\n')
fprintf('1. DW: dwflby2 = %f\n',dwflby2);

% Region 2

dby2by2 = 1;

dbz2by2 = bz1/by1;

dvy2by2 = dvy1by2 + (bx/(wfl*wfl))*(by2 - by1)*dwflby2 ...
	  - (bx/wfl)*(dby2by2 - dby1by2);
	  
dvz2by2 = dvz1by2 + (bx/(wfl*wfl))*(bz2 - bz1)*dwflby2 ...
	  - (bx/wfl)*(dbz2by2 - dbz1by2);
	  
dV2by2 = ((bx/(wfl*by2*by2)) - (by1*V1/(by2*by2)))*dby2by2  ... 
	  + (by1/by2)*dV1by2 + (V1/by2)*dby1by2 ... 
	 + (bx/(by2*wfl*wsl))*(vy2 - vy1)*dwflby2 ...
	 - (bx/(by2*wfl))*(dvy2by2 - dvy1by2);	 
	 
dvx2by2 = dvx1by2 - (V2 - V1)*dwflby2 - wfl*(dV2by2 - dV1by2);

dpg2by2 = dpg1by2 - (by2*dby2by2 - bp1*dbp1by2) ...
	  + (vx2 - vx1)*dwflby2 ...
	  + wfl*(dvx2by2 - dvx1by2);

% compute Alfven wave speed derivative w.r.t. by(2) between regions 2 and 3	  	  
Cal = sqrt(bx*bx/V2);
Cal = -Cal;
dwrdlby2 = -Cal/(2*V2)*dV2by2;
fprintf('2. DW: dwrdlby2 = %f\n',dwrdlby2);	  


% Region 3
dbp3bp2 = 1;

dby3by2 = cos(psi(3))*by2/bperp(2);

dbz3by2 = sin(psi(3))*by2/bperp(2);

dV3by2 = dV2by2;

dvx3by2 = dvx2by2;

dpg3by2 = dpg2by2;

dvy3by2 = dvy2by2 + (bx/(wrdl*wrdl))*(by3 - by2)*dwrdlby2 ...
	  - (bx/wrdl)*(dby3by2 - dby2by2);
	  
dvz3by2 = dvz2by2 + (bx/(wrdl*wrdl))*(bz3 - bz2)*dwrdlby2 ...
	  - (bx/wrdl)*(dbz3by2 - dbz2by2);
	  
% compute slow wave speed derivative w.r.t. by(2) between regions 3 and 4	  
V_d = V3;
pg_d = pg3;
by_d = by3;
bz_d = bz3;
bperp_d = bperp(3);
dV_d = dV3by2;
dpg_d = dpg3by2;
dby_d = dby3by2;

by_u = by4;
dby_u = 1;

bperp_u = bp4;
dbp4bp2 = 0;
dbperp_u = dbp4bp2;
dbperp_d = 1;

wfs = wsl;
	  
wave_speed_derivative;
dwslby2 = dws_dw;
fprintf('3. DW: dwslby2 = %f\n',dwslby2);	  

% Region 4
	  
