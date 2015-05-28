%******************************************************************************* 
%* Program: eigen_mhd.m
%* Description: Computes the eigenvalues and left/right eigenmatrices for
%*              the linearized MHD equations.
%* Author: Andrew Kercher
%* Acknowledgements: This funciton is based on the implementation in
%*                   Athena: an astrophysical MHD code.  See
%*                   https://trac.princeton.edu/Athena/
%* References: 
%*         [1] Dai, W. and Woodward, P. R., "An Approximate Riemann Solver
%*             for Ideal Magnetohydrodynamics", J. Comp. Phys.,
%*             111:354-373, 1994.
%*         [2] Ryu, D. and Jones, T. W., "Numerical Magnetohydrodynamics in
%*             Astrophysics: Algorithm and Tests for One-Dimensional Flow",
%*             ApJ, 442:228-258, March 1995.
%*         [3] J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon
%*             "Athena: A new code for astrophysical MHD", ApJS, (2008),
%*             Appendix B.
%-------------------------------------------------------------------------------
%* Input: neigen = number of eigenvalues (7 for MHD) 
%*         gamma = gas constant
%*          qavg = Averaged state variables (d,vx,vy,vz,h,by,bz,bx)
%*         bperp = Averaged perpendicular magnetic field
%* Output: lambda = eigenvalues
%*          ematL = left eigenmatrix      
%*          ematR = right eigenmatrix      
%******************************************************************************* 

function [lambda,ematL,ematR] = eigen_mhd(neigen,gamma,qavg,bperp)
  
  % initialize eigenvalues and eigenmatrices
  lambda = zeros(neigen,1);
  ematL  = zeros(neigen,neigen);
  ematR  = zeros(neigen,neigen);  
  
  d  = qavg(1);
  vx = qavg(2);
  vy = qavg(3);
  vz = qavg(4); 
  pt = qavg(5);
  bx = qavg(8);
  by = qavg(6);
  bz = qavg(7);
  
%*   bperp = sqrt(by*by + bz*bz);  
  pg = pt - 0.5*(bx*bx + bperp^2);
  v2 = (vx*vx + vy*vy + vz*vz);
  ke = 0.5*v2;
  
  % wave speeds in Lagragian mass coordinates
  [c02,ca2,cperp2,cs2,cf2] = Lagrangian_wave_speeds_mhd(gamma,d,pg,bperp,bx);
  
  % divide by density squared to transform from mass coordinates
  c02 = c02/(d*d);
  ca2 = ca2/(d*d);
  cperp2 = cperp2/(d*d);  
  cs2 = cs2/(d*d);
  cf2 = cf2/(d*d);  
  
  c0 = sqrt(c02);  
  ca = sqrt(ca2);
  cperp = sqrt(cperp2);  
  cs = sqrt(cs2);
  cf = sqrt(cf2);  

  % eigenvalues
  lambda(1) = vx - cf;
  lambda(2) = vx - ca;
  lambda(3) = vx - cs;
  lambda(4) = vx;
  lambda(5) = vx + cs;
  lambda(6) = vx + ca;
  lambda(7) = vx + cf;  
  
  % calulate the sign of bx
  if abs(bx) > 0
    sgnbx = bx/abs(bx);
  else
    sgnbx = 1;    
  end
  
  sgnbp = 1;
  if by > 0
    sgnbp = 1;
  elseif (by == 0) && (bz > 0)
    sgnbp = 1;
  elseif by < 0
    sgnbp = -1;
  elseif (by == 0) && (bz < 0)
    sgnbp = -1;
  end
  
  % define paramaters for computation of left/right eigenmatrix
  if (bperp == 0) && (c02 == ca2)
    alpha_f = 1;
    alpha_s = 1;  
  else
    alpha_f = sqrt(cf2 - ca2)/sqrt(cf2 - cs2);
    alpha_s = sqrt(cf2 - c02)/sqrt(cf2 - cs2);  
  end

  if bperp > 0
    beta_y = by/bperp;
    beta_z = bz/bperp;  
  else
    beta_y = 1/sqrt(2);
    beta_z = 1/sqrt(2);    
  end
  
  gfr = alpha_f*cf2/(gamma - 1) + alpha_f*cf*vx ... 
        - alpha_s*ca*sgnbx*(beta_y*vy + beta_z*vz) ...
        + (gamma - 2)*alpha_f*(cf2 - c02)/(gamma - 1);

  gfl = alpha_f*cf2/(gamma - 1) - alpha_f*cf*vx ... 
        + alpha_s*ca*sgnbx*(beta_y*vy + beta_z*vz) ...
        + (gamma - 2)*alpha_f*(cf2 - c02)/(gamma - 1);

  gsr = alpha_s*cs2/(gamma - 1) + alpha_s*cs*vx ... 
        + alpha_f*c0*sgnbx*(beta_y*vy + beta_z*vz) ...
        + (gamma - 2)*alpha_s*(cs2 - c02)/(gamma - 1);

  gsl = alpha_s*cs2/(gamma - 1) - alpha_s*cs*vx ... 
        - alpha_f*c0*sgnbx*(beta_y*vy + beta_z*vz) ...
        + (gamma - 2)*alpha_s*(cs2 - c02)/(gamma - 1);

  theta_1 = alpha_f^2*c02*(cf2 + (2 - gamma)*c02/(gamma - 1)) ...
            + alpha_s^2*cf2*(cs2 + (2 - gamma)*c02/(gamma - 1));
  
  theta_2 = alpha_f^2*cf*c0*sgnbx + alpha_s^2*cs*ca*sgnbx; 
  
  Qfast_y = alpha_f*beta_y*c0*sgnbx;  
  Qfast_z = alpha_f*beta_z*c0*sgnbx;  

  Qslow_y = alpha_s*beta_y*ca*sgnbx;
  Qslow_z = alpha_s*beta_z*ca*sgnbx;  
  
  % right eigenmatrix
  ematR(1,1) = alpha_f;
  ematR(1,2) = alpha_f*(vx - cf);
  ematR(1,3) = alpha_f*vy + Qslow_y;
  ematR(1,4) = alpha_f*vz + Qslow_z;
  ematR(1,5) = alpha_f*ke + gfl;
  ematR(1,6) = alpha_s*beta_y*cf/sqrt(d);
  ematR(1,7) = alpha_s*beta_z*cf/sqrt(d);

  ematR(2,1) = 0;
  ematR(2,2) = 0;
  ematR(2,3) = beta_z*sgnbx;
  ematR(2,4) = -beta_y*sgnbx;
  ematR(2,5) = (vy*beta_z - vz*beta_y)*sgnbx;
  ematR(2,6) = beta_z/sqrt(d);
  ematR(2,7) = -beta_y/sqrt(d);

  ematR(3,1) = alpha_s;
  ematR(3,2) = alpha_s*(vx - cs);
  ematR(3,3) = alpha_s*vy - Qfast_y;
  ematR(3,4) = alpha_s*vz - Qfast_z;
  ematR(3,5) = alpha_s*ke + gsl;
  ematR(3,6) = -alpha_f*beta_y*c02/(sqrt(d)*cf);
  ematR(3,7) = -alpha_f*beta_z*c02/(sqrt(d)*cf);

  ematR(4,1) = 1.0;
  ematR(4,2) = vx;
  ematR(4,3) = vy;
  ematR(4,4) = vz;
  ematR(4,5) = ke;
  ematR(4,6) = 0;
  ematR(4,7) = 0;

  ematR(5,1) = alpha_s;
  ematR(5,2) = alpha_s*(vx + cs);
  ematR(5,3) = alpha_s*vy + Qfast_y;
  ematR(5,4) = alpha_s*vz + Qfast_z;
  ematR(5,5) = alpha_s*ke + gsr;
  ematR(5,6) = ematR(3,6);
  ematR(5,7) = ematR(3,7);
  
  ematR(6,1) = 0;
  ematR(6,2) = 0;
  ematR(6,3) = -ematR(2,3);
  ematR(6,4) = -ematR(2,4);
  ematR(6,5) = -ematR(2,5);
  ematR(6,6) = ematR(2,6);
  ematR(6,7) = ematR(2,7);;
  
  ematR(7,1) = alpha_f;
  ematR(7,2) = alpha_f*(vx + cf);
  ematR(7,3) = alpha_f*vy - Qslow_y;
  ematR(7,4) = alpha_f*vz - Qslow_z;
  ematR(7,5) = alpha_f*ke + gfr;
  ematR(7,6) = ematR(1,6);
  ematR(7,7) = ematR(1,7);  

  % left eigenmatrix
  ematL(1,1) = alpha_f*c02*v2/(4*theta_1) ...
               + (alpha_f*c0*vx*sgnbx - alpha_s*cs*(vy*beta_y) + vz*beta_z)/(2*theta_2);
  ematL(2,1) = -alpha_f*c02*vx/(2*theta_1) - alpha_f*c0*sgnbx/(2*theta_2);
  ematL(3,1) = -alpha_f*c02*vy/(2*theta_1) + alpha_s*beta_y*cs/(2*theta_2);
  ematL(4,1) = -alpha_f*c02*vz/(2*theta_1) + alpha_s*beta_z*cs/(2*theta_2);
  ematL(5,1) =  alpha_f*c02/(2*theta_1);
  ematL(6,1) = sqrt(d)*alpha_s*beta_y*cf*(cs2 + (2 - gamma)*c02/(gamma - 1))/(2*theta_1);
  ematL(7,1) = sqrt(d)*alpha_s*beta_z*cf*(cs2 + (2 - gamma)*c02/(gamma - 1))/(2*theta_1);

  ematL(1,2) = -0.5*sgnbx*(vy*beta_z - vz*beta_y);
  ematL(2,2) =  0;
  ematL(3,2) =  0.5*sgnbx*beta_z;
  ematL(4,2) = -0.5*sgnbx*beta_y;
  ematL(5,2) =  0;
  ematL(6,2) =  0.5*sqrt(d)*beta_z;
  ematL(7,2) = -0.5*sqrt(d)*beta_y;

  ematL(1,3) =  alpha_s*cf2*v2/(4*theta_1) ...
               + (alpha_s*ca*vx*sgnbx - alpha_f*cf*(vy*beta_y) + vz*beta_z)/(2*theta_2);
  ematL(2,3) = -alpha_s*cf2*vx/(2*theta_1) - alpha_s*ca*sgnbx/(2*theta_2);
  ematL(3,3) = -alpha_s*cf2*vy/(2*theta_1) - alpha_f*beta_y*cf/(2*theta_2);
  ematL(4,3) = -alpha_s*cf2*vz/(2*theta_1) - alpha_f*beta_z*cf/(2*theta_2);
  ematL(5,3) =  alpha_s*cf2/(2*theta_1);
  ematL(6,3) = -sqrt(d)*alpha_f*beta_y*cf*(cf2 + (2 - gamma)*c02/(gamma - 1))/(2*theta_1);
  ematL(7,3) = -sqrt(d)*alpha_f*beta_z*cf*(cf2 + (2 - gamma)*c02/(gamma - 1))/(2*theta_1);

  ematL(1,4) =  1 - v2*(alpha_f^2*c02 + alpha_s^2*cf2)/(2*theta_1);
  ematL(2,4) =  vx*(alpha_f^2*c02 + alpha_s^2*cf2)/theta_1;
  ematL(3,4) =  vy*(alpha_f^2*c02 + alpha_s^2*cf2)/theta_1;
  ematL(4,4) =  vz*(alpha_f^2*c02 + alpha_s^2*cf2)/theta_1;
  ematL(5,4) = -(alpha_f^2*c02 + alpha_s^2*cf2)/theta_1;
  ematL(6,4) =  sqrt(d)*alpha_f*alpha_s*beta_y*cf*(cf2 - cs2)/theta_1;
  ematL(7,4) =  sqrt(d)*alpha_f*alpha_s*beta_z*cf*(cf2 - cs2)/theta_1;

  ematL(1,5) =  alpha_s*cf2*v2/(4*theta_1) ...
               - (alpha_s*ca*vx*sgnbx - alpha_f*cf*(vy*beta_y) + vz*beta_z)/(2*theta_2);
  ematL(2,5) = -alpha_s*cf2*vx/(2*theta_1) + alpha_s*ca*sgnbx/(2*theta_2);
  ematL(3,5) = -alpha_s*cf2*vy/(2*theta_1) + alpha_f*beta_y*cf/(2*theta_2);
  ematL(4,5) = -alpha_s*cf2*vz/(2*theta_1) + alpha_f*beta_z*cf/(2*theta_2);
  ematL(5,5) =  ematL(5,3);
  ematL(6,5) =  ematL(6,3);
  ematL(7,5) =  ematL(7,3);

  ematL(1,6) = -ematL(1,2);
  ematL(2,6) =  0;
  ematL(3,6) = -ematL(3,2);
  ematL(4,6) = -ematL(4,2);
  ematL(5,6) =  0;
  ematL(6,6) =  ematL(6,2);
  ematL(7,6) =  ematL(7,2);

  ematL(1,7) = alpha_f*c02*v2/(4*theta_1) ...
               - (alpha_f*c0*vx*sgnbx - alpha_s*cs*(vy*beta_y) + vz*beta_z)/(2*theta_2);
  ematL(2,7) = -alpha_f*c02*vx/(2*theta_1) + alpha_f*c0*sgnbx/(2*theta_2);
  ematL(3,7) = -alpha_f*c02*vy/(2*theta_1) - alpha_s*beta_y*cs/(2*theta_2);
  ematL(4,7) = -alpha_f*c02*vz/(2*theta_1) - alpha_s*beta_z*cs/(2*theta_2);
  ematL(5,7) =  ematL(5,1);
  ematL(6,7) =  ematL(6,1);
  ematL(7,7) =  ematL(7,1);
  
  % force functions to be continuous
  if (c02 > ca2)
    ematR(3,:) = sgnbp*ematR(3,:);
    ematR(5,:) = sgnbp*ematR(5,:);
    ematL(:,3) = sgnbp*ematL(:,3);
    ematL(:,5) = sgnbp*ematL(:,5);
  else
    ematR(1,:) = sgnbp*ematR(1,:);
    ematR(7,:) = sgnbp*ematR(7,:);
    ematL(:,1) = sgnbp*ematL(:,1);
    ematL(:,7) = sgnbp*ematL(:,7);
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
