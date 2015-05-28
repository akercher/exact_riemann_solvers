%******************************************************************************* 
%* Program: Lagrangian_wave_speed_derivative.m
%* Description: Computes derivatives of various wave speeds in mass coordinates.  
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
%*            V = 1/Density.
%*           pg = Gas pressure.
%*        bperp = Perpendicular magnetic field.
%*           bx = x component of  magnetic field.
%*          C02 = Speed of sound squared.
%*          Ca2 = Parallel Alfven speed squared.
%*       Cperp2 = Perpendicular Alfven speed squared.
%* Output: dC02 = derivative Speed of sound squared.
%*         dCa2 = derivative Parallel Alfven speed squared.
%*      dCperp2 = derivative Perpendicular Alfven speed squared.
%*         dCs2 = derivative Slow wave speed squared.
%*         dCf2 = derivative Fast wave speed squared.
%******************************************************************************* 

function [dC02,dCa2,dCperp2,dCs2,dCf2] = Lagrangian_wave_speed_derivative(gamma,V,pg,...
						  bperp,bx,dV,dpg,dbperp,C02,Ca2,Cperp2)

  % calculate derivatives of speeds squared
  dC02 = (gamma/V)*dpg - (gamma*pg/(V*V))*dV;

  dCa2 = -(bx*bx/(V*V))*dV;
  
  dCperp2 = (2*bperp/V)*dbperp - (bperp*bperp/(V*V))*dV;

  radical = sqrt((C02 + Ca2 + Cperp2)^2 - 4*C02*Ca2);      
  dradical = (C02 + Ca2 + Cperp2)*(dC02 + dCa2 + dCperp2) ...
             - 2*C02*dCa2 - 2*Ca2*dC02;
  dradical = dradical/radical;

  dCf2 = 0.5*(dC02 + dCa2 + dCperp2 + dradical);
  dCs2 = 0.5*(dC02 + dCa2 + dCperp2 - dradical);
