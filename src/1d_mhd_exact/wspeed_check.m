%******************************************************************************* 
%* Program: wspeed_check.m
%* Description: Checks to make sure initial guess satisfies shock
%*              conditions.  We need wfl < wal < wsl < wsr < war < wfr. 
%* Author: Andrew Kercher 
%* References: 
%*         [1] Dai, W. and Woodward, P. R., "An Approximate Riemann Solver
%*             for Ideal Magnetohydrodynamics", J. Comp. Phys.,
%*             111:354-373, 1994.
%*         [2] Ryu, D. and Jones, T. W., "Numerical Magnetohydrodynamics in
%*             Astrophysics: Algorithm and Tests for One-Dimensional Flow",
%*             ApJ, 442:228-258, March 1995.
%-------------------------------------------------------------------------------
%******************************************************************************* 

%----------------------------------------------------------------------
% 1a. Left fast shocks/rarefactions 
%----------------------------------------------------------------------  
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
fs_shock_speed_mhd; 
if wfast2 < 0
  error('ShockSpeed:Imfastl','left fast shock speed has imaginary part...');
end
wfl = -wfast;				% shock speed

% calculate state downstream of left going fast wave
if pstates(5,2) < pstates(5,1) && 1 == 1
  %----------------------------------------------------------------------
  % fast rarefaction (left)
  %----------------------------------------------------------------------  
  pstates(:,2) = fs_rarefaction(gamma,-1,psi(1),bperp(1),bperp(2),...
                                pstates(:,1),bx,'fast');
else
  %----------------------------------------------------------------------
  % shock jump conditions
  %----------------------------------------------------------------------  
  pstates(:,2) = fs_jump(gamma,wfl,psi(1),bperp(2),pstates(:,1),bx);  

end
states(:,2) = prim2cons_mhd_1d(nunks,1,gamma,pstates(:,2),bx);    

%----------------------------------------------------------------------
% 1b. Right fast shocks/rarefactions 
%----------------------------------------------------------------------  
% upstream variables wave traveling to the right
d_u = pstates(1,end);
pg_u = pstates(5,end);  
by_u = pstates(6,end);  
bz_u = pstates(7,end);    
bperp_u = bperp(end);

% downstream variables for fast wave traveling to the right  
bperp_d = bperp(7);			% for RJ calc.
by_d = by7;				% for DW calc.

% compute fast shock speed for wave traveling to the right  
fs_shock_speed_mhd;
if wfast2 < 0
  error('ShockSpeed:Imfastr','right fast shock speed has imaginary part...');
end

wfr = wfast;				% shock speed

% calculate state downstream of right going fast wave
if pstates(5,end-1) < pstates(5,end) && 1 == 1
  %----------------------------------------------------------------------
  % fast rarefaction (right)
  %----------------------------------------------------------------------  
  pstates(:,end-1) = fs_rarefaction(gamma,1,psi(end),bperp(end),...
                                    bperp(end-1),pstates(:,end),bx,'fast');
else
  %----------------------------------------------------------------------
  % shock jump conditions
  %----------------------------------------------------------------------  
  pstates(:,end-1) = fs_jump(gamma,wfr,psi(end),bperp(end-1),pstates(:,end),bx);
end

states(:,end-1) = prim2cons_mhd_1d(nunks,1,gamma,pstates(:,end-1),bx);  

%----------------------------------------------------------------------
% 2a. Left rotaional discontinuity: wave speed is equal to the Alfven 
%       speed in Lagragian coordinates.
%----------------------------------------------------------------------  
% calculate state downstream of left going rotational discontinuity
wal = sqrt(pstates(1,2)*bx*bx);
wal = -wal;  

% check left Alfven speed
if wal < wfl
  % decrease bperp(2)
  bperp(2) = bperp(2)/2;
  return;
end  

pstates(:,3) = rd_jump(wal,psi(3),pstates(:,2),bx);
states(:,3) = prim2cons_mhd_1d(nunks,1,gamma,pstates(:,3),bx);

%----------------------------------------------------------------------
% 2b. Right rotaional discontinuity: wave speed is equal to the Alfven 
%       speed in Lagragian coordinates.
%----------------------------------------------------------------------  
% calculate state downstream of right going rotational discontinuity
war = sqrt(pstates(1,7)*bx*bx);

% check right Alfven speed
if wfr < war
  % decrease bperp(7)
  bperp(7) = bperp(7)/2;
  return;
end  

pstates(:,6) = rd_jump(war,psi(6),pstates(:,7),bx);
states(:,6) = prim2cons_mhd_1d(nunks,1,gamma,pstates(:,6),bx);

%----------------------------------------------------------------------
% 3a. Left slow shocks/rarefactions
%----------------------------------------------------------------------  
% upstream variables slow wave traveling to the left
d_u = pstates(1,3);
pg_u = pstates(5,3);  
by_u = pstates(6,3);  
bz_u = pstates(7,3);    
bperp_u = bperp(3);

% downstream variables for slow wave traveling to the left  
bperp_d = bperp(4);			% for RJ calc.
by_d = by4;				% for DW calc.

fs_shock_speed_mhd;
if wslow2 < 0
  error('ShockSpeed:Imslowl','left slow shock speed has imaginary part...');
end
wsl = -wslow;				% shock speed

% check left slow speed
if wsl <= wfl
  % decrease bperp(4)
  bperp(4) = bperp(4)/4;
  return;
end  

if wsl <= wal
  % decrease bperp(3)
  bperp(3) = bperp(3)/4;
  return;
end  


% calculate state downstream of left going slow wave
if pstates(5,4) < pstates(5,3) && 1 == 1
  %----------------------------------------------------------------------
  % slow rarefaction (left)
  %----------------------------------------------------------------------  
  pstates(:,4) = fs_rarefaction(gamma,-1,psi(3),bperp(3),bperp(4),...
                                pstates(:,3),bx,'slow');
else
  %----------------------------------------------------------------------
  % shock jump conditions
  %----------------------------------------------------------------------
  pstates(:,4) = fs_jump(gamma,wsl,psi(3),bperp(4),pstates(:,3),bx);
end

states(:,4) = prim2cons_mhd_1d(nunks,1,gamma,pstates(:,4),bx);


%----------------------------------------------------------------------
% 3b. Right slow shocks/rarefactions
%----------------------------------------------------------------------  
% upstream variables slow wave traveling to the right
d_u = pstates(1,6);
pg_u = pstates(5,6);  
by_u = pstates(6,6);  
bz_u = pstates(7,6);    
bperp_u = bperp(6);

% downstream variables for slow wave traveling to the left
% across contact discontinuity bperp(4) = bperp(5)
by5 = by4;

bperp_d = bperp(5);			% for RJ calc.
by_d = by5;				% for DW calc.

% compute slow shock speed for wave traveling to the right
fs_shock_speed_mhd;
if wslow2 < 0
  error('ShockSpeed:Imslowr','right slow shock speed has imaginary part...');
end

wsr = wslow;				% shock speed

% check right slow speed
if wfr <= wsr
  % decrease bperp(4)
  bperp(4) = bperp(4)/2;
  return;
end  

if war <= wsr
  % decrease bperp(3)
  bperp(7) = bperp(7)/4;
  return;
end  

% calculate state downstream of right going slow wave
if pstates(5,5) < pstates(5,6) && 1 == 1
  %----------------------------------------------------------------------
  % slow rarefaction (right)
  %----------------------------------------------------------------------  
  pstates(:,5) = fs_rarefaction(gamma,1,psi(6),bperp(6),bperp(5),...
                                pstates(:,6),bx,'slow');
else
  %----------------------------------------------------------------------
  % shock jump conditions
  %----------------------------------------------------------------------  
  pstates(:,5) = fs_jump(gamma,wsr,psi(6),bperp(5),pstates(:,6),bx);

end

states(:,5) = prim2cons_mhd_1d(nunks,1,gamma,pstates(:,5),bx);  

wcond = 1;

