%******************************************************************************* 
%* Program: fs_rarefaction.m
%* Description: Computes the downstream state variables (qdown) using a set
%*              of differential equations derived from the MHD
%*              characteristic equations [3].
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
%*                (d0,vx0,vy0,vz0,pg0,by0,bz0)
%*            bx = x component of  magnetic field
%*          type = type of wave: slow or fast
%* Output: qdown = primitive variables downstream 
%*                (d1,vx1,vy1,vz1,pg1,by1,bz1)
%******************************************************************************* 

function qdown = fs_rarefaction(gamma,dir,psi,bperp0,bperp1,qup,bx,type)

  % upstream variables
  d0 = qup(1);
  vx0 = qup(2);
  vy0 = qup(3);
  vz0 = qup(4);
  pg0 = qup(5);
  by0 = qup(6);
  bz0 = qup(7);
  pt0 = pg0 + 0.5*(bx*bx + by0*by0 + bz0*bz0);

  % 4-stage Runge-Kutta scheme
  h = bperp1 - bperp0;
  y_0 = qup;
  y_0(1) = sqrt(gamma*d0*pg0);

  % first step  
  y_k = y_0;  
  frk4 = ode_wave_func(gamma,dir,psi,bx,y_k,type);    
  rk1 = h*frk4;
  
  % second step 
  y_k = y_0 + 0.5*rk1;
  frk4 = ode_wave_func(gamma,dir,psi,bx,y_k,type);  
  rk2 = h*frk4;
  
  % third step 
  y_k = y_0 + 0.5*rk2;
  frk4 = ode_wave_func(gamma,dir,psi,bx,y_k,type);    
  rk3 = h*frk4;
  
  % forth step 
  y_k = y_0 + rk3;
  frk4 = ode_wave_func(gamma,dir,psi,bx,y_k,type);    
  rk4 = h*frk4;
  
  y_k = y_0 + (rk1 + 2*rk2 + 2*rk3 + rk4)/6;

%*   fprintf('[%d] rk1 = %f rk2 = %f rk3 = %f rk4 = %f\n',0,rk1(1),rk2(1),rk3(1),rk4(1));
%*   fprintf('h = %f y_0 = %f y_k = %f \n',h, y_0(1),y_k(1));
  
  qdown = y_k;  
  d1 = y_k(1)*y_k(1)/(gamma*y_k(5));
  qdown(1) = d1;

  