%******************************************************************************* 
%* Program: Lagrangian_wave_speeds_mhd.m
%* Description: Computes characteristic speeds for plane symmetric MHD
%*              equations in mass coordinates.  
%* Author: Andrew Kercher 
%* References: 
%*         [1] Dai, W. and Woodward, P. R., "An Approximate Riemann Solver
%*             for Ideal Magnetohydrodynamics", J. Comp. Phys.,
%*             111:354-373, 1994.
%*         [2] Ryu, D. and Jones, T. W., "Numerical Magnetohydrodynamics in
%*             Astrophysics: Algorithm and Tests for One-Dimensional Flow",
%*             ApJ, 442:228-258, March 1995.
%-------------------------------------------------------------------------------
%* Input: gamma = Gas constant.
%*            d = Density.
%*           pg = Gas pressure.
%*        bperp = Perpendicular magnetic field.
%*           bx = x component of  magnetic field.
%* Output:  C02 = Speed of sound squared.
%*          Ca2 = Parallel Alfven speed squared.
%*       Cperp2 = Perpendicular Alfven speed squared.
%*          Cs2 = Slow wave speed squared.
%*          Cf2 = Fast wave speed squared.
%******************************************************************************* 

function [C02,Ca2,Cperp2,Cs2,Cf2] = Lagrangian_wave_speeds_mhd(gamma,d,pg,bperp,bx)
  
  %calculate speeds
  C02 = gamma*d*pg;
  Ca2 = d*bx*bx;
  Cperp2 = d*bperp*bperp;

  Cf2_Cs2 = sqrt((C02 + Ca2 + Cperp2)^2 - 4*C02*Ca2);
  Cf2 = 0.5*(C02 + Ca2 + Cperp2 + Cf2_Cs2);
  Cs2 = 0.5*(C02 + Ca2 + Cperp2 - Cf2_Cs2);  
  