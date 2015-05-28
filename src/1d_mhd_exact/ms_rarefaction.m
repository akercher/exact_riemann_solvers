%******************************************************************************* 
%* Program: ms_rarefaction.m
%* Description: Computes the downstream state variables (qdown) using a set
%*              of differential equations derived from the MHD
%*              characteristic equations.
%* Author: Andrew Kercher 
%* References: 
%*         [1] Dai, W. and Woodward, P. R., "An Approximate Riemann Solver
%*             for Ideal Magnetohydrodynamics", J. Comp. Phys.,
%*             111:354-373, 1994.
%*         [2] Ryu, D. and Jones, T. W., "Numerical Magnetohydrodynamics in
%*             Astrophysics: Algorithm and Tests for One-Dimensional Flow",
%*             ApJ, 442:228-258, March 1995.
%*         [3] Jeffrey, A., "Magnetohydrodynamics", Oliver & Boyd, 1966. 
%-------------------------------------------------------------------------------
%* The downstream (post-shock state) is determined from the upstream
%* (pre-shock) state and the downstream perpendicular magnetic field (bperp1). 
%* Input:  gamma = gas constant.
%*           dir = wave direction (left: -1, right: 1)
%*           psi = rotation angle 
%*        bperp0 = perpendicular magnetic field upstream region 
%*        bperp1 = perpendicular magnetic field downstream region 
%*           qup = primitive variables upstream 
%*                 (d0,vx0,vy0,vz0,pg0,by0,bz0)
%*            bx = x component of  magnetic field
%* Output: qdown = primitive variables downstream 
%*                (d1,vx1,vy1,vz1,pg1,by1,bz1)
%******************************************************************************* 

function qdown = ms_rarefaction(gamma,dir,psi,bperp0,bperp1,qup,bx)

  % downstream variables
  d0 = qup(1);
  vx0 = qup(2);
  vy0 = qup(3);
  vz0 = qup(4);
  pg0 = qup(5);
  by0 = qup(6);
  bz0 = qup(7);
  pt0 = pg0 + 0.5*(bx*bx + by0*by0 + bz0*bz0);
  
  % calculate states upstream of shocks
  
  % transverse magnetic field components
  by1 = cos(psi)*bperp1;
  bz1 = sin(psi)*bperp1;
  
  % Riemann invariants and constants
  R0 = pg0/d0^gamma;
  R1 = bperp0/d0;
  K1 = sqrt(gamma*R0);
  K2 = (R1/K1)^2;

  % compute variables 
  d1 = bperp1/R1;			% eq. 3.65 [3]
  pg1 = pg0*(d1/d0)^gamma;
  vx1 = vx0 + dir*2*(K1/K2)*((1 + K2*d1^(1/3))^(3/2) - (1 + K2*d0^(1/3))^(3/2));
  
  qdown(1) = d1;
  qdown(2) = vx1;
  qdown(3) = 0;
  qdown(4) = 0;
  qdown(5) = pg1;
  qdown(6) = by1;  
  qdown(7) = bz1;    
  
