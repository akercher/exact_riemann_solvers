%******************************************************************************* 
%* Program: ms_rstates_mhd.m
%* Description: Computes MHD Riemann states from initial left/right
%*              states. 
%* Author: Andrew Kercher 
%* References: 
%*         [1] Dai, W. and Woodward, P. R., "An Approximate Riemann Solver
%*             for Ideal Magnetohydrodynamics", J. Comp. Phys.,
%*             111:354-373, 1994.
%*         [2] Ryu, D. and Jones, T. W., "Numerical Magnetohydrodynamics in
%*             Astrophysics: Algorithm and Tests for One-Dimensional Flow",
%*             ApJ, 442:228-258, March 1995.
%*         [3] Toro, E. F., "Riemann Solvers and Numerical Methods for
%*             Fluid Dynamics", Springer-Verlag, 2nd ed., 1999.
%-------------------------------------------------------------------------------
%*
%* Case: abs(bx) = 0.  States in regions 4 and 5 must satisfy jump
%*   conditions for a tangential discontinuity, that is only the densitiey
%*   and total pressure should have discontinuity.  This leads to the
%*   following system
%*             vx4(bperp2,bperp4,psi) = vx5(bperp4,bperp7,psi)
%*                        pt4(bperp2) = pt5(bperp7)
%*   that can be solved iteratively for bperp2,bperp7.  The
%*   iterativel procedure is a modified version of that given by
%*   equations (11a) -(11d) of [1] where by has been replaced with bperp.  
%*
%******************************************************************************* 

%----------------------------------------------------------------------
% magnetosonic shocks: RJ method
%----------------------------------------------------------------------  
  
% calculate state left of tangential discontinuity (region 2)
% upstream variables fast wave traveling to the left
d_u = pstates(1,1);
pg_u = pstates(5,1);  
by_u = pstates(6,1);  
bz_u = pstates(7,1);    
bperp_u = bperp(1);

% downstream variables for fast wave traveling to the left  
bperp_d = bperp(2);			% for RJ calc.
by_d = by2;				% for DW calc.

% compute fast shock speed for wave traveling to the left
% fast wave speed equals magnetosonic wave speed
fs_shock_speed_mhd;
wmsl = -wfast;				% shock speed
%* Mal = abs(wmsl)/Cf;			% Mach number  
Cfl = Cf;				% wave speed  
Cfavgl = Cf;

if pstates(5,2) <= pstates(5,1)
  %----------------------------------------------------------------------
  % magnetosonic fast rarefaction (left)
  %----------------------------------------------------------------------  
  pstates(:,2) = ms_rarefaction(gamma,-1,psi(1),bperp(1),bperp(2),pstates(:,1),bx);
  states(:,2) = prim2cons_mhd_1d(nunks,1,gamma,pstates(:,2),bx);  
else
  %----------------------------------------------------------------------
  % magnetosonic shock jump conditions
  %----------------------------------------------------------------------  
  pstates(:,2) = fs_jump(gamma,wmsl,psi(1),bperp(2),pstates(:,1),bx);  
  states(:,2) = prim2cons_mhd_1d(nunks,1,gamma,pstates(:,2),bx);  
end

% upstream variables wave traveling to the right
d_u = pstates(1,end);
pg_u = pstates(5,end);  
by_u = pstates(6,end);  
bz_u = pstates(7,end);    
bperp_u = bperp(end);  

% downstream variables for fast wave traveling to the right  
bperp_d = bperp(end-1);			% for RJ calc.
by_d = by7;				% for DW calc.

% compute fast shock speed for wave traveling to the left
% fast wave speed equals magnetosonic wave speed
fs_shock_speed_mhd;  
wmsr = wfast;				% shock speed
%* Mar = abs(wmsr)/Cf;			% Mach number  
Cfr = Cf;				% wave speed  
Cfavgr = Cf;

% calculate state right of tangential discontinuity (region 3)
if pstates(5,end-1) <= pstates(5,end)
  %----------------------------------------------------------------------
  % magnetosonic fast rarefaction (right)
  %----------------------------------------------------------------------        
  pstates(:,end-1) = ms_rarefaction(gamma,1,psi(end),bperp(end),bperp(end-1),...
                                    pstates(:,end),bx);
  states(:,end-1) = prim2cons_mhd_1d(nunks,1,gamma,pstates(:,3),bx);
else
  %----------------------------------------------------------------------
  % magnetosonic shock jump conditions
  %----------------------------------------------------------------------  
  pstates(:,end-1) = fs_jump(gamma,wmsr,psi(end),bperp(end-1),pstates(:,end),bx);  
  states(:,end-1) = prim2cons_mhd_1d(nunks,1,gamma,pstates(:,3),bx);  
end

% set states to left of TD
pstates(:,3) = pstates(:,2);
pstates(:,4) = pstates(:,2);
states(:,3) = states(:,2);
states(:,4) = states(:,2);
%* pstates(:,3:4) = pstates(:,2);
%* states(:,3:4) = states(:,2);

% set states to right of TD
pstates(:,5) = pstates(:,end-1);
pstates(:,6) = pstates(:,end-1);
states(:,5) = states(:,end-1);
states(:,6) = states(:,end-1);
%* pstates(:,5:6) = pstates(:,end-1);
%* states(:,5:6) = states(:,end-1);

