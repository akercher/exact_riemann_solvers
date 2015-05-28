%******************************************************************************* 
%* Program: ms_jacob_p.m
%* Description: Computes Jacobian (derivative matrix) for magnetosonic
%*              shocks using total pressure formulation given by DW [1].  
%* Author: Andrew Kercher 
%* References: 
%*         [1] Dai, W. and Woodward, P. R., "An Approximate Riemann Solver
%*             for Ideal Magnetohydrodynamics", J. Comp. Phys.,
%*             111:354-373, 1994.
%-------------------------------------------------------------------------------
%* Input:   wms = magnetosonic shock speed
%*        qdown = primative variables downstream
%*                (d0,vx0,vy0,vz0,pg0,by0,bz0)
%*          qup = primative variables upstream
%*                (d1,vx1,vy1,vz1,pg1,by1,bz1)
%* Output:  DJ  = derivative matrix for Newton iteration
%******************************************************************************* 

function DJ = ms_jacob_p(gamma,wms,qdown,qup)

  % downstream variables
  d0 = qdown(1);
  pg0 = qdown(5);  
  by0 = qdown(6);  
  bz0 = qdown(7);    
  
  % total pressures 
  pt0 = pg0 + 0.5*(by0^2 + bz0^2);
  pt1 = qup(5) + 0.5*(qup(6)^2 + qup(7)^2);  

  % difference of total pressure
  pt_d = pt1 - pt0;  
  
% speeds in Lagrangian mass coordiantes  
  C02 = gamma*pg0*d0;
  Cperp2 = d0*(by0*by0 + bz0*bz0);
  Cms2 = C02 + Cperp2;  
  
  % dw/dp
  G1 = 0.5*(gamma + 3)*d0*pt_d;
  G2 = 2*d0*d0*((gamma + 1)*pt_d + 2*gamma*pt0)*pt_d;
  dG1 = 0.5*(gamma + 3)*d0;
  dG2 = 2*d0*d0*(gamma + 1)*pt_d + G2/pt_d;
  
  dCms2 = 0;%gamma*d0;
  
  radical = sqrt((Cms2 + G1)^2 - G2);
  dradical =  (2*(Cms2 + G1)*dG1 - dG2)/(2*radical);
  
  dwms2 = 0.5*(dCms2 + dG1 + dradical);
  
  % dwms/dp
  dwdp = dwms2/(2*wms);

  % dvx/dw
  dvxdw = (pt_d)/(wms*wms);
  
  % dvx/dp
  dvxdp = 1/wms; 
  
  DJ = dvxdp + dvxdw*dwdp; 
  