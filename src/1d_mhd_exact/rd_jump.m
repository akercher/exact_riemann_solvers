%******************************************************************************* 
%* Program: rd_jump.m
%* Description: Computes the upstream state variables (qup) with the MHD
%*              jump conditions for rotational discontinuities.
%* Author: Andrew Kercher 
%* References: 
%*         [1] Dai, W. and Woodward, P. R., "An Approximate Riemann Solver
%*             for Ideal Magnetohydrodynamics", J. Comp. Phys.,
%*             111:354-373, 1994.
%*         [2] Ryu, D. and Jones, T. W., "Numerical Magnetohydrodynamics in
%*             Astrophysics: Algorithm and Tests for One-Dimensional Flow",
%*             ApJ, 442:228-258, March 1995.
%-------------------------------------------------------------------------------
%* The upstream (post-shock state) is determined from the downstream
%* (pre-shock) state and the upstream rotation angle (psi). 
%* Input:    ws = fast/slow shock speed
%*          psi = rotation angle
%*           q0 = primitive variables downstream 
%*                (d0,vx0,vy0,vz0,pg0,by0,bz0)
%*           bx = x component of  magnetic field
%* Output:  q1 = primitivevariables upstream 
%*                (d1,vx1,vy1,vz1,pg1,by1,bz1)
%******************************************************************************* 

function q1 = rd_jump(ws,psi,q0,bx)

  % downstream variables
  V0 = 1/q0(1);
  vx0 = q0(2);
  vy0 = q0(3);
  vz0 = q0(4);
  pg0 = q0(5);
  by0 = q0(6);
  bz0 = q0(7);
  bperp0 = sqrt(by0*by0 + bz0*bz0);
  
  % calculate states upstream of shocks
  % rotational discontinuities are noncompressive
  V1 = V0;
  vx1 = vx0;
  pg1 = pg0;

  % transverse magnetic field components  
  by1 = cos(psi)*bperp0;
  bz1 = sin(psi)*bperp0;

  vy1 = vy0 - bx*(by1 - by0)/ws;  
  vz1 = vz0 - bx*(bz1 - bz0)/ws;
  
%*   pt1 = pg1 + 0.5*(bx*bx + by1*by1 + bz1*bz1);    
%*   en1 = en0*V0/V1 + (vx1*pt1 - vx0*pt0)/(V1*ws) - bx*(bx*(vx1 - vx0) ... 
%*           + vy1*by1 - vy0*by0 + vz1*bz1 - vz0*bz0)/(V1*ws);    
  
  q1(1) = 1/V1;
  q1(2) = vx1;
  q1(3) = vy1;
  q1(4) = vz1;
  q1(5) = pg1;
  q1(6) = by1;  
  q1(7) = bz1;    
  
  
  