%******************************************************************************* 
%* Program: rarefaction_waves.m
%* Description: Computes value for rarefaction waves at sampled points
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
%******************************************************************************* 

%---------------------------------------------------------------  
% left fast rarefaction wave.
%---------------------------------------------------------------  
if ~isempty(lfr_cells) && length(lfr_cells) >= 1
  
  dir = -1;
  type = 'fast';

  bperp0 = bperp(1);		% upstream 
  bperp1 = bperp(2);		% downstream
    
  psi_k = psi(1);
  
  ns = length(lfr_cells);
  ds = (bperp1-bperp0)/(ns+1);

  s = (bperp0:ds:bperp1)';
  
  for is=1:ns;

    % 4-stage Runge-Kutta scheme    
    icell = lfr_cells(is);
    h = ds;
    s_0 = s(is+1); 
    if icell == 1  
      q_0 = pstates(:,1);
    else
      q_0 = pstate(:,icell-1);
    end
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

    pstate(:,icell) = q_k;    
    d1 = q_k(1)*q_k(1)/(gamma*q_k(5));
    pstate(1,icell) = d1;

  end
end  

%---------------------------------------------------------------  
% left slow rarefaction wave.
%---------------------------------------------------------------  
if ~isempty(lsr_cells) && length(lsr_cells) >= 1
  
  dir = -1;
  type = 'slow';

  bperp0 = bperp(end/2-1);		% upstream 
  bperp1 = bperp(end/2);		% downstream
    
  psi_k = psi(end/2);
  
  ns = length(lsr_cells);
  ds = (bperp1-bperp0)/(ns+1);

  s = (bperp0:ds:bperp1)';
  
  for is=1:ns;

    % 4-stage Runge-Kutta scheme    
    icell = lsr_cells(is);
    h = ds;
    s_0 = s(is+1);  
    q_0 = pstate(:,icell-1);
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

    pstate(:,icell) = q_k;    
    d1 = q_k(1)*q_k(1)/(gamma*q_k(5));
    pstate(1,icell) = d1;
%*     keyboard
  end
end

%---------------------------------------------------------------  
% right slow rarefaction wave.
%---------------------------------------------------------------  
if ~isempty(rsr_cells) && length(rsr_cells) >= 1
  
  dir = 1;
  type = 'slow';

  bperp0 = bperp(end/2+2);		% upstream 
  bperp1 = bperp(end/2+1);		% downstream
    
  psi_k = psi(end/2+1);
  
  ns = length(rsr_cells);
  ds = (bperp1-bperp0)/(ns+1);

  s = (bperp0:ds:bperp1)';
  
  for is=1:ns;

    % 4-stage Runge-Kutta scheme    
    icell = rsr_cells((ns+1) - is);
    h = ds;
    s_0 = s(is+1);  
    q_0 = pstate(:,icell+1);
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

    pstate(:,icell) = q_k;    
    d1 = q_k(1)*q_k(1)/(gamma*q_k(5));
    pstate(1,icell) = d1;

  end
end

%---------------------------------------------------------------  
% right fast rarefaction wave.
%---------------------------------------------------------------  
if ~isempty(rfr_cells) && length(rfr_cells) >= 1

  dir = 1;
  type = 'fast';
  
  bperp0 = bperp(end);			% upstream 
  bperp1 = bperp(end-1);		% downstream
    
  psi_k = psi(end);
  
  ns = length(rfr_cells);
  ds = (bperp1-bperp0)/(ns+1);

  s = (bperp0:ds:bperp1)';
  
  for is=1:ns;

    % 4-stage Runge-Kutta scheme    
    icell = rfr_cells((ns+1) - is);
    h = ds;
    s_0 = s(is+1);  
    q_0 = pstate(:,icell+1);
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

    pstate(:,icell) = q_k;    
    d1 = q_k(1)*q_k(1)/(gamma*q_k(5));
    pstate(1,icell) = d1;
  
  end
end
