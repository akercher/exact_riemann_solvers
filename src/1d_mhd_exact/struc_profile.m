%******************************************************************************* 
%* Program: struc_profile.m
%* Description: Computes profile for various structures.
%* Author: Andrew Kercher 
%* References: 
%-------------------------------------------------------------------------------
%******************************************************************************* 

%---------------------------------------------------------------  
% postitions of diffent regions
%---------------------------------------------------------------
xr = zeros(nregio,2);
xr(1,1) = grid_1d.xmin;
if grid_1d.xpos(1,1) < grid_1d.xmin
  xr(1,1) = grid_1d.xpos(1,1);
end
xr(1,2) = grid_1d.xpos(1,1);
for i=2:nregio-1
  xr(i,1) = grid_1d.xpos(i-1,2);
  xr(i,2) = grid_1d.xpos(i,1);  
end  
xr(end,1) = grid_1d.xpos(end,2);
xr(end,2) = grid_1d.xmax;
if grid_1d.xpos(end,2) > grid_1d.xmax
  xr(end,2) = grid_1d.xpos(end,2);
end

%---------------------------------------------------------------  
% left fast structure
%---------------------------------------------------------------  
if abs(grid_1d.xpos(1,1) - grid_1d.xpos(1,2)) < tol
  % shock
  ns = 2;
  xfl = [grid_1d.xpos(1,1),grid_1d.xpos(1,2)];
  qfl = zeros(nunks,ns);
  qfl(:,1) = pstates(:,1);
  qfl(:,2) = pstates(:,2);  
else
  
  % left fast rarefaction  
  dir = -1;
  type = 'fast';

  bperp0 = bperp(1);		% upstream 
  bperp1 = bperp(2);		% downstream
    
  psi_k = psi(1);
  
  smin = 0;
  smax = -log(pstates(1,2)/pstates(1,1));
  
  ns = nrp;
  ds = (bperp1-bperp0)/(ns-1);
  s = (bperp0:ds:bperp1)';  
  
  dxs = abs(grid_1d.xpos(1,2) - grid_1d.xpos(1,1))/(ns-1);
  xfl = grid_1d.xpos(1,1):dxs:grid_1d.xpos(1,2);
  qfl = zeros(nunks,ns);
  qfl(:,1) = pstates(:,1);
  qfl(:,end) = pstates(:,2);  
  for is=2:ns-1;

    % 4-stage Runge-Kutta scheme    
    s_0 = s(is); 
    q_0 = qfl(:,is-1);

    d0 = q_0(1);
    pg0 = q_0(5);
    q_0(1) = sqrt(gamma*d0*pg0);

    % first step  
    q_k = q_0;  

    frk4 = ode_wave_func(gamma,dir,psi_k,bx,q_k,type);  
    rk1 = ds*frk4;

    % second step 
    q_k = q_0 + 0.5*rk1;
    frk4 = ode_wave_func(gamma,dir,psi_k,bx,q_k,type);  
    rk2 = ds*frk4;

    % third step 
    q_k = q_0 + 0.5*rk2;
    frk4 = ode_wave_func(gamma,dir,psi_k,bx,q_k,type);  
    rk3 = ds*frk4;

    % forth step 
    q_k = q_0 + rk3;
    frk4 = ode_wave_func(gamma,dir,psi_k,bx,q_k,type);  
    rk4 = ds*frk4;

    q_k = q_0 + (rk1 + 2*rk2 + 2*rk3 + rk4)/6;
  
    qfl(:,is) = q_k;    
    d1 = q_k(1)*q_k(1)/(gamma*q_k(5));
    qfl(1,is) = d1;

    % get rid of signed zero
    if abs(qfl(7,is)) < tol 
      qfl(7,is) = 0;
    end
  
  end  
end
ufl = zeros(nunks,ns);  
ufl = prim2cons_mhd_1d(nunks,ns,gamma,qfl,bx);  


%---------------------------------------------------------------  
% left Alfven structure
%---------------------------------------------------------------  
ns = 2;
xal = [grid_1d.xpos(2,1),grid_1d.xpos(2,2)];
qal = zeros(nunks,ns);
qal(:,1) = pstates(:,2);
qal(:,2) = pstates(:,3);
ual = zeros(nunks,ns);
ual = prim2cons_mhd_1d(nunks,ns,gamma,qal,bx);

%---------------------------------------------------------------  
% left slow structure
%---------------------------------------------------------------  
if abs(grid_1d.xpos(3,1) - grid_1d.xpos(3,2)) < tol
  % shock
  ns = 2;
  xsl = [grid_1d.xpos(3,1),grid_1d.xpos(3,2)];
  qsl = zeros(nunks,ns);
  qsl(:,1) = pstates(:,3);
  qsl(:,2) = pstates(:,4);
else
  % left slow rarefaction  
  dir = -1;
  type = 'slow';

  bperp0 = bperp(end/2-1);		% upstream 
  bperp1 = bperp(end/2);		% downstream
    
  psi_k = psi(end/2);
  
  ns = nrp;
  ds = (bperp1-bperp0)/(ns-1);
  s = (bperp0:ds:bperp1)';
  
  dxs = abs(grid_1d.xpos(3,2) - grid_1d.xpos(3,1))/(ns-1);
  xsl = grid_1d.xpos(3,1):dxs:grid_1d.xpos(3,2);
  qsl = zeros(nunks,ns);
  qsl(:,1) = pstates(:,3);
  qsl(:,end) = pstates(:,4);  
  for is=2:ns-1;

    % 4-stage Runge-Kutta scheme    
    s_0 = s(is); 
    q_0 = qsl(:,is-1);

    d0 = q_0(1);
    pg0 = q_0(5);
    q_0(1) = sqrt(gamma*d0*pg0);

    % first step  
    q_k = q_0;  

    frk4 = ode_wave_func(gamma,dir,psi_k,bx,q_k,type);  

    rk1 = ds*frk4;

    % second step 
    q_k = q_0 + 0.5*rk1;
    frk4 = ode_wave_func(gamma,dir,psi_k,bx,q_k,type);  
    rk2 = ds*frk4;

    % third step 
    q_k = q_0 + 0.5*rk2;
    frk4 = ode_wave_func(gamma,dir,psi_k,bx,q_k,type);  

    rk3 = ds*frk4;

    % forth step 
    q_k = q_0 + rk3;
    frk4 = ode_wave_func(gamma,dir,psi_k,bx,q_k,type);  
    rk4 = ds*frk4;

    q_k = q_0 + (rk1 + 2*rk2 + 2*rk3 + rk4)/6;

    qsl(:,is) = q_k;    
    d1 = q_k(1)*q_k(1)/(gamma*q_k(5));
    qsl(1,is) = d1;

    % get rid of signed zero
    if abs(qsl(7,is)) < tol 
      qsl(7,is) = 0;
    end
  
  end % left slow rarefaction  
end
usl = zeros(nunks,ns);
usl = prim2cons_mhd_1d(nunks,ns,gamma,qsl,bx);  

%---------------------------------------------------------------  
% Contact Discontinuity
%---------------------------------------------------------------  
ns = 2;
xcd = [grid_1d.xpos(4,1),grid_1d.xpos(4,2)];
qcd = zeros(nunks,ns);
qcd(:,1) = pstates(:,4);
qcd(:,2) = pstates(:,5);
ucd = zeros(nunks,ns);
ucd = prim2cons_mhd_1d(nunks,ns,gamma,qcd,bx);  

%---------------------------------------------------------------  
% right slow structure
%---------------------------------------------------------------  
if abs(grid_1d.xpos(5,1) - grid_1d.xpos(5,2)) < tol
  % shock
  ns = 2;
  xsr = [grid_1d.xpos(5,1),grid_1d.xpos(5,2)];
  qsr = zeros(nunks,ns);
  qsr(:,1) = pstates(:,5);
  qsr(:,2) = pstates(:,6);
else
  % right slow rarefaction  
  dir = 1;
  type = 'slow';

  bperp0 = bperp(end/2+2);		% upstream 
  bperp1 = bperp(end/2+1);		% downstream
    
  psi_k = psi(end/2+1);
  
  ns = nrp;
  ds = (bperp1-bperp0)/(ns-1);
  s = (bperp0:ds:bperp1)';
  
  dxs = abs(grid_1d.xpos(5,2) - grid_1d.xpos(5,1))/(ns-1);
  xsr = grid_1d.xpos(5,1):dxs:grid_1d.xpos(5,2);
  qsr = zeros(nunks,ns);
  qsr(:,1) = pstates(:,5);
  qsr(:,end) = pstates(:,6);
  
  for is=2:ns-1;

    ix = ns - is + 1;
  
    % 4-stage Runge-Kutta scheme    
    s_0 = s(is); 
    q_0 = qsr(:,ix+1);

    d0 = q_0(1);
    pg0 = q_0(5);
    q_0(1) = sqrt(gamma*d0*pg0);

    % first step  
    q_k = q_0;  

    frk4 = ode_wave_func(gamma,dir,psi_k,bx,q_k,type);  

    rk1 = ds*frk4;

    % second step 
    q_k = q_0 + 0.5*rk1;
    frk4 = ode_wave_func(gamma,dir,psi_k,bx,q_k,type);  
    rk2 = ds*frk4;

    % third step 
    q_k = q_0 + 0.5*rk2;
    frk4 = ode_wave_func(gamma,dir,psi_k,bx,q_k,type);  

    rk3 = ds*frk4;

    % forth step 
    q_k = q_0 + rk3;
    frk4 = ode_wave_func(gamma,dir,psi_k,bx,q_k,type);  
    rk4 = ds*frk4;

    q_k = q_0 + (rk1 + 2*rk2 + 2*rk3 + rk4)/6;

    qsr(:,ix) = q_k;    
    d1 = q_k(1)*q_k(1)/(gamma*q_k(5));
    qsr(1,ix) = d1;

    % get rid of signed zero
    if abs(qsr(7,ix)) < tol 
      qsr(7,ix) = 0;
    end
  
  end % right slow rarefaction
end
usr = zeros(nunks,ns);
usr = prim2cons_mhd_1d(nunks,ns,gamma,qsr,bx);  

%---------------------------------------------------------------  
% right Alfven structure
%---------------------------------------------------------------  
xar = [grid_1d.xpos(6,1),grid_1d.xpos(6,2)];
ns = 2;
qar = zeros(nunks,ns);
qar(:,1) = pstates(:,6);
qar(:,2) = pstates(:,7);
uar = zeros(nunks,ns);
uar = prim2cons_mhd_1d(nunks,ns,gamma,qar,bx);  

%---------------------------------------------------------------  
% right fast structure
%---------------------------------------------------------------  
if abs(grid_1d.xpos(7,1) - grid_1d.xpos(7,2)) < tol
  % shock
  ns = 2;
  xfr = [grid_1d.xpos(7,1),grid_1d.xpos(7,2)];
  qfr = zeros(nunks,ns);
  qfr(:,1) = pstates(:,7);
  qfr(:,2) = pstates(:,8);  
else
  % right fast rarefaction  
  dir = 1;
  type = 'fast';

  bperp0 = bperp(end);		% upstream 
  bperp1 = bperp(end-1);		% downstream
    
  psi_k = psi(end);
  
  ns = nrp;
  ds = (bperp1-bperp0)/(ns-1);
  s = (bperp0:ds:bperp1)';
  
  dxs = abs(grid_1d.xpos(7,2) - grid_1d.xpos(7,1))/(ns-1);
  xfr = grid_1d.xpos(7,1):dxs:grid_1d.xpos(7,2);
  qfr = zeros(nunks,ns);
  qfr(:,1) = pstates(:,7);
  qfr(:,end) = pstates(:,8);
  
  for is=2:ns-1;

    ix = ns - is + 1;

    % 4-stage Runge-Kutta scheme    
    s_0 = s(is); 
    q_0 = qfr(:,ix + 1);

    d0 = q_0(1);
    pg0 = q_0(5);
    q_0(1) = sqrt(gamma*d0*pg0);

    % first step  
    q_k = q_0;  

    frk4 = ode_wave_func(gamma,dir,psi_k,bx,q_k,type);  

    rk1 = ds*frk4;

    % second step 
    q_k = q_0 + 0.5*rk1;
    frk4 = ode_wave_func(gamma,dir,psi_k,bx,q_k,type);  
    rk2 = ds*frk4;

    % third step 
    q_k = q_0 + 0.5*rk2;
    frk4 = ode_wave_func(gamma,dir,psi_k,bx,q_k,type);  

    rk3 = ds*frk4;

    % forth step 
    q_k = q_0 + rk3;
    frk4 = ode_wave_func(gamma,dir,psi_k,bx,q_k,type);  
    rk4 = ds*frk4;

    q_k = q_0 + (rk1 + 2*rk2 + 2*rk3 + rk4)/6;

    qfr(:,ix) = q_k;    
    d1 = q_k(1)*q_k(1)/(gamma*q_k(5));
    qfr(1,ix) = d1;

    % get rid of signed zero
    if abs(qfr(7,ix)) < tol 
      qfr(7,ix) = 0;
    end
  
  end % right fast rarefaction
end
ufr = zeros(nunks,ns);
ufr = prim2cons_mhd_1d(nunks,ns,gamma,qfr,bx);    
