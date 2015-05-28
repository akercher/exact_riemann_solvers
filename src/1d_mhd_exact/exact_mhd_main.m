%******************************************************************************* 
%* Program: exact_mhd_main.m
%* Description: Exact nonlinear Riemann sovler for 1D Magnetohydrodynamic
%*              equations.
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
%*                                t
%*                               /|\ Contact or Tangential Discontinuity
%*                                |
%*          Slow (vx-cs)\         |         /Slow (vx+cs)
%*     Rotation .        \    R4  |  R5    /        . Rotation
%*      (vx-ca)   .   R3  \       |       /  R6   .    (vx+ca)
%*                  .      \      |      /      .
%* Fast.              .     \     |     /     .              .Fast 
%*  (vx-cf) .    R2     .    \    |    /    .     R7    .      (vx+cf)
%*               .        .   \   |   /   .       .
%*                    .     .  \  |  /  .    .
%*           R1            .  . \ | / .  .           R8
%*   /___________________________\|/___________________________\ x
%*   \                                                         /
%*     The 7 possible waves and 8 possible states from a the MHD Riemann
%*     problem [2].  The speeds are given in terms of the characteristic
%*     speeds for the plane symmetric ideal MHD equations.  The slow
%*     speed is cs, the fast speed is cf and the Alfven speed is ca.
%*
%* The MHD jump conditiions in the mass coordinate, dm = rho*dx, are
%* given as (eq. numbers refer to [2])
%*   W*[V]    = -[vx]                                (3.1) 
%*   W*[vx]   = -[P - Bx*Bx]                         (3.2)  
%*   W*[vy]   = -Bx*[By]                             (3.3)
%*   W*[vz]   = -Bx*[Bz]                             (3.4)
%*   W*[V*By] = -Bx*[vy]                             (3.5)
%*   W*[V*Bz] = -Bx*[vz]                             (3.6)
%*   W*[V*E]  = [vx*P] - Bx*[Bx*vx + By*vy + Bz*vz]  (3.7)
%* where P = p + (Bx*Bx +By*By + Bz*Bz)/2 and W = -rho*vx is the
%* Largragian speed of the discontinuity and V = 1/rho [2].  When moving
%* in the -x (or -m in the mass coordinate) direction, the speed W is
%* negative. The braket [Q] denotes the difference between the downstream
%* and upstream states of a quantity Q, i.e [Q] = Qd - Qu.  
%* 
%* Rotational discontinuities, set [vx] = 0 in equations (3.1) and (3.2).    
%*
%* Contact discontinuities, set W = 0 in equations (3.1)-(3.7).
%*
%* Tangential discontinuities, special case of a contact discontinuity
%* where the normal component of the magnetic field is zero (Bn = 0).  In
%* one-dimension this occurs when Bx = 0.  Tangential components of
%* magnetic field and flow velocity are not necessarily equal across the
%* discontinuity.  Only equations (3.1), (3.2) and (3.7), modified for a
%* contact discontinuity (W = 0) need to be satisfied.  
%*
%* Fast and slow shocks have the following relationship between
%* downstream and upstream states (equation numbers refer to [1])
%*   Bzd*[By] = Byd*[Bz]  (4a)
%*   Bzd*[vy] = Byd*[vz]  (5)
%* Thus, from eq. 4a we have
%*   Bzd/Byd = Bzu/Byu    (4b)
%* meaning the orientations of the transverse magnetic field are either
%* exactly the same or opposite.
%-------------------------------------------------------------------------------
%*
%* Computational procedure:
%*
%* 1. Guess initial states in all eight regions.  Region 1 is initial
%*    left state and region 8 is initial right state.
%* 2. Starting with region 1, find the upstream (region 2) state using
%*    the above jump conditions and the tangential magnetic field
%*    components obtained from the inital guess as follows:
%*    a. calculate fast shock speed, eq. (3.8) of [2]
%*    b. Solve eqs 3.3 and 3.4 for upstream values of vy and vz in
%*       region 2, with upstream values of By and Bz obtained from
%*       initial guess. 
%*    b. Use eq. 3.5 or 3.6 to find upstream (region 2) value of V.
%*    c. Use eq. 3.1 to find downstream (region 2) value of vx.
%*    d. Use eq. 3.2 to find downstream (region 2) value of P.
%*    e. Use eq. 3.7 to find downstream (region 2) value of E.
%*    Follow the same procedure to obtain the state in region 7 from the
%     upstream intial state in region 8.
%* 3. Using the rotation angles, defined by tan(phi) = Bz/By, and Alfven
%*    speed from regions 3 and 6, solve the appropriate jump conditions
%*    for rotational discontinuities to find the states in regions 3 and
%*    6 from the upstream states in regions 2 and 7 respectively.
%* 4. Using states from regions 3 and 6, compute state in region 4 and 5
%*    with the appropriate jump conditions and speed for slow shocks.
%*    This is the same as step 2.
%* 5. Check to see if the jump conditions for a contact discontinuities
%*    is satisfied.  If not, improve initial guess and repeat process
%*    starting with step 2. 
%*
%******************************************************************************* 

%-------------------------------------------------------------------------------
% Initialize arrays
%-------------------------------------------------------------------------------
if strcmp(run_type,'new')
  states = zeros(nunks,nregio);		% conservative state variables
  pstates = zeros(nunks,nregio);	% primitive state variables
  bperp = zeros(nregio,1);		% perpendicular magnetic field
  psi = zeros(nregio,1);		% rotation angle
  xpos = zeros(neigen,2);		% position of each structure
  if (lstate == 1) && (rstate == 8)
    grid_1d.xmin = x1min;
    grid_1d.xmax = x1max;
    grid_1d.xd = xd0;
  end
  grid_1d.xpos = zeros(neigen,2);	% position of each structure  
  Lx1 = x1max - x1min;			% length of domain
end

latex_ws = [];				% store wave speeds and residuals

pstates(1,1) = dl;
pstates(2,1) = vxl;
pstates(3,1) = vyl;
pstates(4,1) = vzl;
pstates(5,1) = pgl;
pstates(6,1) = byl;
pstates(7,1) = bzl;

pstates(1,end) = dr;
pstates(2,end) = vxr;
pstates(3,end) = vyr;
pstates(4,end) = vzr;
pstates(5,end) = pgr;
pstates(6,end) = byr;
pstates(7,end) = bzr;

%-------------------------------------------------------------------------------
% Non-dimensionalize initial states
%-------------------------------------------------------------------------------
if ~strcmp(state_dimensions,'none')
  [pstates,states,grid_1d,tf,bx,state_dimensions] = ...
                states_nondim(nunks,nregio,gamma,pstates,grid_1d,tf,bx);  
end

% define perpendicular magnetic field in left (region 1) and right states (region 8)
bperp(1) = sqrt(pstates(6,1)^2 + pstates(7,1)^2);
bperp(end) = sqrt(pstates(6,end)^2 + pstates(7,end)^2);

states(:,1) = prim2cons_mhd_1d(nunks,1,gamma,pstates(:,1),bx);
states(:,end) = prim2cons_mhd_1d(nunks,1,gamma,pstates(:,end),bx);

% set pressure so initial calc. is shock
pstates(5,2:4) = large_n;
pstates(5,5:end-1) = large_n;

% define rotation angle in left (region 1) and right states (region 8)

if ~exist('psil','var')
  psil = atan2(bzl,byl);  
end
if ~exist('psir','var')
  psir = atan2(bzr,byr);  
end

psi(1:4) = psil;
psi(5:end) = psir;  

%-------------------------------------------------------------------------------
% Print problem information
%-------------------------------------------------------------------------------
% finish if calc is for flux computation
% create and open log file, pwd = current directory
fname = sprintf('%s/bin/log.%s.out',pwd,prob);
flog = fopen(fname,'w');

% create and open log file for wave speeds, pwd = current directory
fname = sprintf('%s/bin/log_ws.%s.out',pwd,prob);
flog_ws = fopen(fname,'w');

% print to log file
fprintf(flog,'Problem id        : %s\n',prob);
fprintf(flog,'Jacobian method   : %s\n',Jacobian_method);
fprintf(flog,'Initialization    : %s\n',init_guess);
fprintf(flog,'Error Calc.       : %s\n',error_type);
fprintf(flog,'tol               : %0.2e\n',tol);
fprintf(flog,'Relax. fac.       : %0.2e\n',rfac_k);
fprintf(flog,'gamma             : %0.2f\n',gamma);
fprintf(flog,'Bn                : %0.2f\n',bx);
fprintf(flog,'Num. sample cells : %d\n',ncell);
fprintf(flog,'Final time        : %0.4e\n',tf);

% print to standard output
fprintf('\n');
fprintf(' Problem id        : %s\n',prob);
fprintf(' Jacobian method   : %s\n',Jacobian_method);
fprintf(' Initialization    : %s\n',init_guess);
fprintf(' Error Calc.       : %s\n',error_type);
fprintf(' tol               : %0.2e\n',tol);
fprintf(' Relax. fac.       : %0.2e\n',rfac_k);
fprintf(' gamma             : %0.2f\n',gamma);
fprintf(' Bn                : %0.2f\n',bx);
fprintf(' Num. sample cells : %d\n',ncell);
fprintf(' Final time        : %0.4e\n',tf);

%-------------------------------------------------------------------------------
% Calculate initial guess
%-------------------------------------------------------------------------------

guess_mhd;    

if abs(bx) > 0 
  % update rotation angle and bperp arrays
  bperp(3) = bperp(2);
  bperp(5) = bperp(4);
  bperp(6) = bperp(7);  
  psi(4:6) = psi(3);
  
else
  bperp(3:4) = bperp(2);
  bperp(5:6) = bperp(7);
  psi(2:4) = psi(1);
  psi(5:7) = psi(end);  
end

% check for switch-structures
if abs(bperp(1)) < tol
  psi(1:2) = psi(3);
end

if abs(bperp(end)) < tol  
  psi(end-1:end) = psi(3);  
end

%-------------------------------------------------------------------------------
% Compute the solution, here is where the majic happens.
%-------------------------------------------------------------------------------
% solve for Riemann states in all regions
if strcmp(Jacobian_method,'nested')
  % use nested derivatives to compute Jacobian [1].
  rsolver_mhd;
elseif strcmp(Jacobian_method,'approx')
  % approx. Jacobian with finite differences.
  rsolver_mhd_fd;
else
  fprintf(' Error: invalid method for Jacobian computation.\n');
  fprintf('        Valid options: nested or approx.\n');
  return;
end

if ~strcmp(run_type,'flux_calc') 
fprintf(' Num. iterations   : %d\n',kiter);
fprintf(' Error             : %0.2e\n',err_k);
fprintf('\n');
end

% set to zero if less than tol.
for i=1:nregio
  if abs(pstates(2,i)) < tol
    pstates(2,i) = 0;
  end
  if abs(pstates(3,i)) < tol
    pstates(3,i) = 0;
  end
  if abs(pstates(4,i)) < tol
    pstates(4,i) = 0;
  end
  if abs(pstates(6,i)) < tol
    pstates(6,i) = 0;
  end
  if abs(pstates(7,i)) < tol
    pstates(7,i) = 0;
  end
end  

entropy = pstates(5,:)./(pstates(1,:).^(gamma));

% update conservative variables
states = prim2cons_mhd_1d(nunks,nregio,gamma,pstates,bx);

% wave speeds
wspd = [wfl;wal;wsl;wsr;war;wfr];

% save to file
save('./bin/vars','bperp','psi','by2','by4','by7','rfac_k','pstates','state_dimensions');

if strcmp(sampling,'on')
%-------------------------------------------------------------------------------
% Sampling
%-------------------------------------------------------------------------------
% grid variables
Lx = grid_1d.xmax - grid_1d.xmin;
xd0 = grid_1d.xd;
dx = (Lx)/ncell;
xcell = (0.5:ncell-0.5)'/ncell;
xcell = xcell*Lx + grid_1d.xmin;
pstate = zeros(nunks,ncell);		% final primitive state profile
state = zeros(nunks,ncell);		% final conservative state profile
psi_x = zeros(ncell,1);			% rotation angle

lfr_cells = [];			% fast left fast rarefaction
lsr_cells = [];			% left slow rarefaction
rsr_cells = [];			% right slow rarefaction
rfr_cells = [];			% right fast rarefaction

% calculate solution at time tf
for icell=1:ncell
  
  % find sampled wave speed
  wsi = (xcell(icell) - xd0)/tf;

  % compoute solution (di,vxi,vyi,vzi,pgi,byi,bzi) at point (x,t)
  sample_mhd;

end

% calculate rarefaction waves for sampled data
rarefaction_waves;

for icell = 1:ncell
  % get rid of signed zero
  if abs(pstate(7,icell)) < tol 
    pstate(7,icell) = 0;
  end
  % compute rotation angle in degrees  
  psi_x(icell) = atan2(pstate(7,icell),pstate(6,icell))*180/pi;    
end  

% compute conservative variables
state = prim2cons_mhd_1d(nunks,ncell,gamma,pstate,bx);

end

%-------------------------------------------------------------------------------
% Create profile of structures and regions
%-------------------------------------------------------------------------------
% calulate the position of various structures
if ~exist('nrp','var')
  nrp = 256;
end
[grid_1d.xpos,vspd] = struc_pos(gamma,pstates,bperp,wspd,grid_1d.xd,tf,bx);
struc_profile;

%-------------------------------------------------------------------------------
% Convert to quatities with dimensions
%-------------------------------------------------------------------------------
%* states_dim;

% call plotting routine
rhoL = pstates(1,1);
vxL = pstates(2,1);
vyL = pstates(3,1);
vzL = pstates(4,1);
pL = pstates(5,1);
ByL = pstates(6,1);
BzL = pstates(7,1);
BxL = bx;

rhoR = pstates(1,end);
vxR = pstates(2,end);
vyR = pstates(3,end);
vzR = pstates(4,end);
pR = pstates(5,end);
ByR = pstates(6,end);
BzR = pstates(7,end);
BxR = bx;

% Short format: export data to text file, pwd = current directory
fname = sprintf('%s/dat/%s.out',pwd,prob);
fid = fopen(fname,'w');
if abs(bx) > 0
fprintf(fid,'wave speeds: -fast, -Alfven, -slow, slow, Alfven, fast\n');
fprintf(fid,'% -0.4e % -0.4e % -0.4e % -0.4e % -0.4e % -0.4e\n',wfl,wal,wsl,wsr,war,wfr);
else
fprintf(fid,'wave speeds: -Magnetosonic, Magnetosonic\n');
fprintf(fid,'% -0.4e % -0.4e\n',wmsl,wmsr);  
end
fprintf(fid,'\n');
fprintf(fid,'state variables: rho, vx, vy, vz, by, bz, pg\n');
for i=1:nregio
  fprintf(fid,'% -0.4e % -0.4e % -0.4e % -0.4e % -0.4e % -0.4e % -0.4e\n',pstates(1,i),pstates(2,i),pstates(3,i),pstates(4,i),pstates(6,i),pstates(7,i),pstates(5,i));  
end
fprintf(fid,'\n');
fprintf(fid,'structure positions: xstart, xend\n');
for i=1:neigen
  fprintf(fid,'% -0.4e % -0.4e\n',xpos(i,1),xpos(i,2));
end

% close file
fclose(fid);

% Long format: export data to text file, pwd = current directory
fname = sprintf('%s/dat/%s_long.out',pwd,prob);
fid = fopen(fname,'w');
if abs(bx) > 0
fprintf(fid,'wave speeds: -fast, -Alfven, -slow, slow, Alfven, fast\n');  
fprintf(fid,'% -0.*e % -0.*e % -0.*e % -0.*e % -0.*e % -0.*e\n',fmt_pre,wfl,fmt_pre,wal,fmt_pre,wsl,fmt_pre,wsr,fmt_pre,war,fmt_pre,wfr);
else
fprintf(fid,'wave speeds: -Magnetosonic, Magnetosonic\n');  
fprintf(fid,'% -0.*e % -0.*e\n',fmt_pre,wmsl,fmt_pre,wmsr);  
end
fprintf(fid,'\n');
fprintf(fid,'state variables: rho, vx, vy, vz, by, bz, pg\n');
for i=1:nregio
  fprintf(fid,'% -0.*e % -0.*e % -0.*e % -0.*e % -0.*e % -0.*e % -0.*e\n',fmt_pre,pstates(1,i),fmt_pre,pstates(2,i),fmt_pre,pstates(3,i),fmt_pre,pstates(4,i),fmt_pre,pstates(6,i),fmt_pre,pstates(7,i),fmt_pre,pstates(5,i));  
end
fprintf(fid,'\n');
fprintf(fid,'structure positions: xstart, xend\n');
for i=1:neigen
  fprintf(fid,'% -0.*e % -0.*e\n',fmt_pre,xpos(i,1),fmt_pre,xpos(i,2));
end

% close files
%* fclose(fid);
fclose(flog);
fclose(flog_ws);
