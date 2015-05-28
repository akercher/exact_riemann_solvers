%******************************************************************************* 
%* Program: Lagrangian_shock_speeds_mhd.m
%* Description: Computes fast/slow shock speed in Lagrangian mass coordinates.
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
%* Variables are defined in Lagrangian mass coordinates.
%*  Input:    Cs2 = Slow characteristic speed squared.
%*            Cf2 = Fast characteristic speed squared.
%*             S0 = Wave speed coeffient, eq. (3.9) [2].
%*             S1 = Wave speed coeffient, eq. (3.10) [2].
%*             S2 = Wave speed coeffient, eq. (3.11) [2].
%*  Output: wfast2 = Fast wave speed squared, eq. (3.8) of [2]. 
%*          wslow2 = Slow wave speed squared, eq. (3.8) of [2]. 
%*           errio = return error if negative under radical
%******************************************************************************* 

function [wfast2,wslow2,errio] = Lagrangian_shock_speeds_mhd(Cs2,Cf2,S0,S1,S2)
  
  errio = 0;
  
  % fast/slow shock speeds squared, eq. (9) of [1] and eq. 3.8 [2]
  CsCfS1 = (Cs2 + Cf2 + S1);
  radical = CsCfS1*CsCfS1 - 4*(1+S0)*(Cs2*Cf2 - S2);
  
  if radical < 0
    errio = 1
  end
  
  radical = sqrt(radical);  
  wfast2 = 0.5*(CsCfS1 + radical)/(1+S0);
  wslow2 = 0.5*(CsCfS1 - radical)/(1+S0);
  
  