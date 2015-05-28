%******************************************************************************* 
%* Program: ode_wave_func.m
%* Description: Computes functions for rhs of system of odes for
%*              fast/slow rarefaction waves in ideal MHD (see
%*              eqs. 3.12-3.17 of [2]).
%* Author: Andrew Kercher 
%* References: 
%*         [1] Dai, W. and Woodward, P. R., "An Approximate Riemann Solver
%*             for Ideal Magnetohydrodynamics", J. Comp. Phys.,
%*             111:354-373, 1994.
%*         [2] Ryu, D. and Jones, T. W., "Numerical Magnetohydrodynamics in
%*             Astrophysics: Algorithm and Tests for One-Dimensional Flow",
%*             ApJ, 442:228-258, March 1995.
%*         [3] Jeffrey, A., "Magnetohydrodynamics", Oliver & Boyd, 1966. 
%*         [4] Brio, M. & Wu, C. C., "An upwind difference shceme for the
%              equation of ideal magnetohydrodynamics."  J. Comp. Phys., 
%*             75:400-422, 1988. 
%-------------------------------------------------------------------------------
%* Computes values for rhs of eqs. 3.12-3.17 [2] to be used to solve the
%* system:  d(q_k)/dbperp = fode.  
%* Input:  gamma = gas constant.
%*           dir = wave direction (left: -1, right: 1)
%*           psi = rotation angle 
%*            bx = x component of  magnetic field
%*           q_k = vector of unknowns, primitive variable with Lagranian
%*                 speed of sound replacing density.
%*                 (C0,vx0,vy0,vz0,pg0,by0,bz0)
%*          type = fast or slow wave, character string
%* Output:  fode = rhs functions
%******************************************************************************* 

function fode = ode_wave_func(gamma,dir,psi,bx,q_k,type)
  
  fode = zeros(size(q_k));
  
  C0 = q_k(1);
%*   d = q_k(1);  
  vx = q_k(2);
  vy = q_k(3);
  vz = q_k(4);
  pg = q_k(5);  
  by = q_k(6);
  bz = q_k(7);  
  
  d = C0*C0/(gamma*pg);  
%*   C0 = (gamma*d*pg);    
  sqrtd = sqrt(d);
  bperp = sqrt(by*by + bz*bz);

  [C02,Ca2,Cperp2,Cs2,Cf2] = Lagrangian_wave_speeds_mhd(gamma,d,pg,bperp,bx);
%*   C0 = sqrt(C02);
  Ca = sqrt(Ca2);  
  Cperp = sqrt(Cperp2);
  Cf = sqrt(Cf2);  
  Cs = sqrt(Cs2);    

  if strcmp(type,'fast')
  
    fode(1) = -(gamma+1)*sqrtd*Cperp*Cs2/(2*C0*(Cs2 - Ca2));
%*     fode(1) = d*Cs2*(Cf2 - Ca2)/(C02*Cperp*Ca2);    
    fode(2) = -dir*Cperp*Ca2/(sqrtd*Cf*(Cs2 - Ca2));
    fode(3) = -dir*cos(psi)*Ca/(sqrtd*Cf);
    fode(4) = -dir*sin(psi)*Ca/(sqrtd*Cf);  
    fode(5) = Cs2*(Cf2 - Ca2)/(sqrtd*Cperp*Ca2);  
    fode(6) = cos(psi);  
    fode(7) = sin(psi);    
  
  elseif strcmp(type,'slow')
    fode(1) = -(gamma+1)*sqrtd*Cperp*Cf2/(2*C0*(Cf2 - Ca2));
%*     fode(1) = d*Cf2*(Cs2 - Ca2)/(C02*Cperp*Ca2);    
    fode(2) = -dir*Cperp*Ca2/(sqrtd*Cs*(Cf2 - Ca2));
    fode(3) = -dir*cos(psi)*Ca/(sqrtd*Cs);
    fode(4) = -dir*sin(psi)*Ca/(sqrtd*Cs);  
    fode(5) = Cf2*(Cs2 - Ca2)/(sqrtd*Cperp*Ca2);  
    fode(6) = cos(psi);  
    fode(7) = sin(psi);      
  
%*     printf('frk[0] = %f C0 = %f Ct = %f Cf2 = %f\n',fode(1),C0,Cperp,Cf2);

  else
    fprintf(' Error: wave type must be either fast or slow.\n');
    fprintf(' Function: ode_wave_func.\n');  
    return;
  
  end

  
  