%******************************************************************************* 
%* Program: rsolver_mhd_fd.m
%* Description: Computes MHD Riemann states from initial left/right
%*              states. Approximates the Jacobian with finite differences.
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
%* Case I: abs(bx) > 0.  States in regions 4 and 5 must satisfy jump
%*   conditions for a contact discontinuity, that is only the densities and
%*   energies should be differnet.  This leads to the following system
%*             vx4(bperp2,bperp4,psi) = vx5(bperp4,bperp7,psi)
%*             vy4(bperp2,bperp4,psi) = vy5(bperp4,bperp7,psi)
%*             vz4(bperp2,bperp4,psi) = vz5(bperp4,bperp7,psi)
%*             pg4(bperp2,bperp4,psi) = pg5(bperp4,bperp7,psi)
%*   that can be solved iteratively for bperp2, bperp4,bperp7,psi.  The
%*   iterativel procedure is given by equations (11a) -(11d) of [2] where
%*   by has been replaced with bperp.  
%* Case II: abs(bx) = 0.  States in regions 4 and 5 must satisfy jump
%*   conditions for a tangential discontinuity, that is the normal
%*   velocity and total pressures will be equal.  In this case, the
%*   following equation must be solved.
%*             vx4(pt4,W-(pt4)) = vx5(pt4,W+(pt4)).
%*   The iterative procedure is can be found in [2].
%*
%******************************************************************************* 

% change in iteration variables 
if abs(bx) > 0
  
  % system of nonlinear equations
  fcd = zeros(4,1);			% vector of functions at CD
  
  xold = zeros(4,1);    
  xvec = zeros(4,1);			% vector of independent variables
  xvec(1) = bperp(2);
  xvec(2) = bperp(4);
  xvec(3) = bperp(7);
  xvec(4) = psi(3);  
  
  err = zeros(4,1);
  zavg = zeros(4,1);
  
else
  
  % system of nonlinear equations
  ftd = zeros(2,1);			% vector at TD

  xold = zeros(2,1);    
  xvec = zeros(2,1);			% vector of independent variables  
  xvec(1) = bperp(2);
  xvec(2) = bperp(end-1);
  
  err = zeros(2,1);
  zavg = zeros(2,1);
  
end

if abs(bx) > 0

  %----------------------------------------------------------------------
  % Fast shocks, slow shocks and rotational discontinuities
  %----------------------------------------------------------------------  

  kiter = 0;
  
  % store value of independent variables
  xold(:) = xvec(:);

  % initailize solution assuming all shocks, no rarefactions
  fs_rstates_mhd;

  % store shock/wave speeds
  wfl0 = wfl;				% left fast shock speed
  Cfl0 = Cfl;				% left fast wave speed  
  wal0 = wal;				% left rotational discontinuity speed
  wsl0 = wsl;				% left slow shock speed
  Csl0 = Csl;				% left slow wave speed  
  wsr0 = wsr;				% right slow shock speed
  Csr0 = Csr;				% right slow wave speed  
  war0 = war;				% right rotational discontinuity speed				  
  wfr0 = wfr;				% right fast shock speed
  Cfr0 = Cfr;				% right fast wave speed
  
  % rotated angles in degrees
  dpsil = (psi(3) - psi(2))*180/pi;
  dpsir = (psi(6) - psi(7))*180/pi;  
  
  % difference in state variables, we want this to converge to zero
  fcd(1) = pstates(2,5) - pstates(2,4);
  fcd(2) = pstates(3,5) - pstates(3,4);
  fcd(3) = pstates(4,5) - pstates(4,4);
  fcd(4) = pstates(5,5) - pstates(5,4) + 0.5*(bperp(5)^2 -  bperp(4)^2);
  
  % zone averaged values
  zavg(1) = 0.5*(pstates(2,5) + pstates(2,4));
  zavg(2) = 0.5*(pstates(3,5) + pstates(3,4));
  zavg(3) = 0.5*(pstates(4,5) + pstates(4,4));
  zavg(4) = 0.5*(pstates(5,5) + pstates(5,4) ... 
            + 0.5*(2*bx^2 + bperp(5)^2 + bperp(4)^2));

  zavg(1) = max(zavg(1),0.5*(Cfavgr+Cfavgl));
  zavg(2) = zavg(1); 
  zavg(3) = zavg(1);   

  % compute Jacobian with finite differences
  fs_jacob_fd;

  dxvec = DJ \ fcd;
  
  xvec(:) = xvec(:) - dxvec;  % ****** solve system and subtract *********

  % check if jump conditions are satisfied at contact discontinuity
  if strcmp(error_type,'zone_avg')
    err(:) = abs(fcd(:))./zavg(:);

  elseif strcmp(error_type,'max_change')
    % compute error as maximum change
    err(:) = 2*abs((xvec(:) - xold(:))./(xvec(:) + xold(:)));
  else
    fprintf('\n\n');
    fprintf('  ERROR: invalid choice for error calculation.\n');
    return;
  end

  err_k = max(err(:));
  
  % revert to previous values
  xvec(:) = xold(:);
  
  kiter = 0;  

  % revert back to stored shock/wave speeds
  wfl = wfl0;				% left fast shock speed
  Cfl = Cfl0;				% left fast wave speed  
  wal = wal0;				% left rotational discontinuity speed
  wsl = wsl0;				% left slow shock speed
  Csl = Csl0;				% left slow wave speed  
  wsr = wsr0;				% right slow shock speed
  Csr = Csr0;				% right slow wave speed  
  war = war0;				% right rotational discontinuity speed
  wfr = wfr0;				% right fast shock speed
  Cfr = Cfr0;				% right fast wave speed

  % rotated angles in degrees
  dpsil = (psi(3) - psi(2))*180/pi;
  dpsir = (psi(6) - psi(7))*180/pi;  

  % write to log file
%*   log_output;
  
%*   fprintf('\n');
  fprintf(flog_ws,'                                Waves Speeds and Rotations\n');  
  fprintf(flog_ws,' -------------------------------------------------------------------------------------------\n');  
  fprintf(flog_ws,'   iter       wfl           wsl           wsr           wfr         dpsil         dpsir\n');
  fprintf(flog_ws,' -------------------------------------------------------------------------------------------\n');
  fprintf(flog_ws,'  % 4d   % -0.4e   % -0.4e   % -0.4e   % -0.4e   % -0.4e   % -0.4e\n',kiter,wfl,wsl,wsr,wfr,dpsil,dpsir);   

  while  (tol <= err_k) && (kiter < maxit)  

    kiter = kiter + 1;

    % update independent variables
    bperp(2:3) = xvec(1);
    bperp(4:5) = xvec(2);
    bperp(6:7) = xvec(3);
    psi(3:6) = xvec(4);

    % update rotation angle and bperp arrays
    by2 = cos(psi(2))*bperp(2);
    by4 = cos(psi(4))*bperp(4);
    by7 = cos(psi(7))*bperp(7);

    % solve jump condition or integrate rarefactions
    fs_rstates_mhd;

    % store shock/wave speeds
    wfl0 = wfl;				% left fast shock speed
    Cfl0 = Cfl;				% left fast wave speed  
    wal0 = wal;			% left rotational discontinuity speed
    wsl0 = wsl;				% left slow shock speed
    Csl0 = Csl;				% left slow wave speed  
    wsr0 = wsr;				% right slow shock speed
    Csr0 = Csr;				% right slow wave speed  
    war0 = war;			% right rotational discontinuity speed				
    wfr0 = wfr;				% right fast shock speed
    Cfr0 = Cfr;				% right fast wave speed

    % rotated angles in degrees
    dpsil = (psi(3) - psi(2))*180/pi;
    dpsir = (psi(6) - psi(7))*180/pi;  

    %----------------------------------------------------------------------
    % Actual iteration step
    %----------------------------------------------------------------------  

    % difference in state variables, we want this to converge to zero
    fcd(1) = pstates(2,5) - pstates(2,4);
    fcd(2) = pstates(3,5) - pstates(3,4);
    fcd(3) = pstates(4,5) - pstates(4,4);
    fcd(4) = pstates(5,5) - pstates(5,4) + 0.5*(bperp(5)^2 -  bperp(4)^2);

    % zone averaged values
    zavg(1) = 0.5*(pstates(2,5) + pstates(2,4));
    zavg(2) = 0.5*(pstates(3,5) + pstates(3,4));
    zavg(3) = 0.5*(pstates(4,5) + pstates(4,4));
    zavg(4) = 0.5*(pstates(5,5) + pstates(5,4) ... 
              + 0.5*(2*bx^2 + bperp(5)^2 + bperp(4)^2));

    zavg(1) = max(zavg(1),0.5*(Cfavgr+Cfavgl));
    zavg(2) = zavg(1);
    zavg(3) = zavg(1);

    % compute Jacobian with finite differences
    fs_jacob_fd;  

    dxtemp = DJ \ fcd(:);  

    % calculate optimal step size
    if(~(max(abs(rfac_k * dxtemp)) < max(abs(dxvec)) ))

      % decrease relaxation factor
      while( ( max(abs(dxvec)) < max(abs(rfac_k * dxtemp))))
        
        % till step size decreases
        rfac_k = 0.5*rfac_k;
        if( rfac_k < rfac_min )
          fprintf('Damping too strong. No convergence.\n');
          return; 
        end
      end % step size is valid
    end

    % update variables
    dxvec = dxtemp;
    xvec(:) = xvec(:) - rfac_k*dxvec;	% ****** subtract inverted Jacobian *********

    % check for switch-structures
    if abs(bperp(1)) < tol
      psi(1:2) = psi(3);
    end

    if abs(bperp(end)) < tol  
      psi(end-1:end) = psi(3);  
    end

    % check if jump conditions are satisfied at contact discontinuity
    if strcmp(error_type,'zone_avg')
      err(:) = abs(fcd(:))./zavg(:);

    elseif strcmp(error_type,'max_change')
      % compute error as maximum change
      err(:) = 2*abs((xvec(:) - xold(:))./(xvec(:) + xold(:)));
    else
      fprintf('\n\n');
      fprintf('  ERROR: invalid choice for error calculation.\n');
      return;
    end

    err_k = max(err(:));
%*     err_k = max(rfac_k*abs(dxvec(:)))

    % increase relaxation factor
    rfac_k = min(1,(rfac_k + 0.5*rfac_k));

    xold(:) = xvec(:);

    % revert back to stored shock/wave speeds
    wfl = wfl0;				% left fast shock speed
    Cfl = Cfl0;				% left fast wave speed  
    wal = wal0;			% left rotational discontinuity speed
    wsl = wsl0;				% left slow shock speed
    Csl = Csl0;				% left slow wave speed  
    wsr = wsr0;				% right slow shock speed
    Csr = Csr0;				% right slow wave speed  
    war = war0;			% right rotational discontinuity speed
    wfr = wfr0;				% right fast shock speed
    Cfr = Cfr0;				% right fast wave speed
  
    % rotated angles in degrees
    dpsil = (psi(3) - psi(2))*180/pi;
    dpsir = (psi(6) - psi(7))*180/pi;  

    % write to log file
    log_output;

    % write to standard output
    fprintf(flog_ws,'  % 4d   % -0.4e   % -0.4e   % -0.4e   % -0.4e   % -0.4e   % -0.4e\n',kiter,wfl,wsl,wsr,wfr,dpsil,dpsir);   

  end

  % compute Mach numbers
  Mafl = abs(wfl)/Cfl;			% Mach number  
  Mafr = abs(wfr)/Cfr;			% Mach number  
  Masl = abs(wsl)/Cfl;			% Mach number  
  Masr = abs(wsr)/Cfr;			% Mach number  

  fprintf(flog_ws,'\n\n');
  fprintf(flog_ws,'                       Mach Numbers\n');
  fprintf(flog_ws,' --------------------------------------------------------\n');
  fprintf(flog_ws,'      Mafl          Mafr          Masl          Masr\n');
  fprintf(flog_ws,' --------------------------------------------------------\n');
  fprintf(flog_ws,'  % -0.4e   % -0.4e   % -0.4e   % -0.4e\n',Mafl,Mafr,Masl,Masr);   
%*   fprintf(flog_ws,'\n');


else

  %----------------------------------------------------------------------
  % magnetosonic shocks: RJ method
  %----------------------------------------------------------------------  

  % store value of independent variables
  xold(:) = xvec(:);
  
  % initailize solution assuming all shocks, no rarefactions
  ms_rstates_mhd;  

  % store shock/wave speeds
  wmsl0 = wmsl;				% left fast shock speed
  Cfl0 = Cfl;				% left fast wave speed  
  wmsr0 = wmsr;				% right fast shock speed
  Cfr0 = Cfr;				% right fast wave speed
  
  % difference in state variables, we want this to converge to zero
  ftd(1) = pstates(2,end-1) - pstates(2,2);
  ftd(2) = pstates(5,end-1) - pstates(5,2) + 0.5*(bperp(end-1)^2 -  bperp(2)^2);
  
  % zone averaged values
  zavg(1) = 0.5*(pstates(2,end-1) + pstates(2,2));
  zavg(2) = 0.5*(pstates(5,end-1) + pstates(5,2) ... 
            + 0.5*(2*bx^2 + bperp(end-1)^2 + bperp(2)^2));

  zavg(1) = max(zavg(1),0.5*(Cfavgr+Cfavgl));
  
  % compute Jacobian with finite differences
  ms_jacob_fd;

  dxvec = DJ \ ftd;
  xvec(:) = xvec(:) - dxvec;  % ****** solve system and subtract *********

  % check if jump conditions are satisfied at contact discontinuity
  if strcmp(error_type,'zone_avg')
    err(:) = abs(ftd(:))./zavg(:);

  elseif strcmp(error_type,'max_change')
    % compute error as maximum change
    err(:) = 2*abs((xvec(:) - xold(:))./(xvec(:) + xold(:)));
  else
    fprintf('\n\n');
    fprintf('  ERROR: invalid choice for error calculation.\n');
    return;
  end

  err_k = max(err(:));

  % revert to previous values
  xvec(:) = xold(:);
  
  kiter = 0;  

  % revert back to stored shock/wave speeds
  wmsl = wmsl0;				% left fast shock speed
  Cfl = Cfl0;				% left fast wave speed  
  wmsr = wmsr0;				% right fast shock speed
  Cfr = Cfr0;				% right fast wave speed

  % write to log file
  log_output;

  % print initail states
%*   print_states_mhd;
  
%*   fprintf('\n');
  fprintf(flog_ws,'     Magnetosonic Waves Speeds\n');
  fprintf(flog_ws,' ----------------------------------\n');  
  fprintf(flog_ws,'   iter      wmsl          wmsr \n');
  fprintf(flog_ws,' ----------------------------------\n');
  fprintf(flog_ws,'  % 4d   % -0.4e   % -0.4e\n',kiter,wmsl,wmsr);       
  
  % iterate to find root
  while (tol <= err_k) && (kiter < maxit)  

    kiter = kiter + 1;

    % update variables
    bperp(2:4) = xvec(1);
    bperp(5:7) = xvec(2);

    % solve jump condition or integrate rarefactions
    ms_rstates_mhd;

    % store shock/wave speeds
    wmsl0 = wmsl;				% left fast shock speed
    Cfl0 = Cfl;				% left fast wave speed  
    wmsr0 = wmsr;				% right fast shock speed
    Cfr0 = Cfr;				% right fast wave speed

    %----------------------------------------------------------------------
    % Actual iteration step
    %----------------------------------------------------------------------  

    % difference in state variables, we want this to converge to zero
    ftd(1) = pstates(2,end-1) - pstates(2,2);
    ftd(2) = pstates(5,end-1) - pstates(5,2) + 0.5*(bperp(end-1)^2 -  bperp(2)^2);
  
    % zone averaged values
    zavg(1) = 0.5*(pstates(2,end-1) + pstates(2,2));
    zavg(2) = 0.5*(pstates(5,end-1) + pstates(5,2) ... 
              + 0.5*(2*bx^2 + bperp(end-1)^2 + bperp(2)^2));

    zavg(1) = max(zavg(1),0.5*(Cfavgr+Cfavgl));
  
    % compute Jacobian with finite differences
    ms_jacob_fd;

    dxtemp = DJ \ ftd(:);  

    % calculate optimal step size
    if(~(max(abs(rfac_k * dxtemp)) < max(abs(dxvec)) ))

      % decrease relaxation factor
      while( ( max(abs(dxvec)) < max(abs(rfac_k * dxtemp))))
        
        % till step size decreases
        rfac_k = 0.5*rfac_k;
        if( rfac_k < rfac_min )
          fprintf('Damping too strong. No convergence.\n');
          return; 
        end
      end % step size is valid
    end

    % update variables
    dxvec = dxtemp;
    xvec(:) = xvec(:) - rfac_k*dxvec;	% ****** subtract inverted Jacobian *********

    % increase relaxation factor
    rfac_k = min(1,(rfac_k + 0.5*rfac_k));

    % check if jump conditions are satisfied at contact discontinuity
    if strcmp(error_type,'zone_avg')
      err(:) = abs(ftd(:))./zavg(:);

    elseif strcmp(error_type,'max_change')
      % compute error as maximum change
      err(:) = 2*abs((xvec(:) - xold(:))./(xvec(:) + xold(:)));
    else
      fprintf('\n\n');
      fprintf('  ERROR: invalid choice for error calculation.\n');
      return;
    end

    err_k = max(err(:));

    xold(:) = xvec(:);

    % revert back to stored shock/wave speeds
    wmsl = wmsl0;				% left fast shock speed
    Cfl = Cfl0;				% left fast wave speed  
    wmsr = wmsr0;				% right fast shock speed
    Cfr = Cfr0;				% right fast wave speed
  
    % write to log file
    log_output;

    % write to standard output
    fprintf(flog_ws,'  % 4d   % -0.4e   % -0.4e\n',kiter,wmsl,wmsr);       

  end
  
%*   fprintf('  %04d   %0.4f   %0.4f\n',kiter,wmsl,wmsr);         
  
  Mal = abs(wmsl)/Cfl;			% Mach number  
  Mar = abs(wmsr)/Cfr;			% Mach number  

  fprintf(flog_ws,'\n\n');  
  fprintf(flog_ws,'         Mach Numbers\n');
  fprintf(flog_ws,' ----------------------------\n');
  fprintf(flog_ws,'       Mal           Mar\n');
  fprintf(flog_ws,' ----------------------------\n');
  fprintf(flog_ws,'  % -0.4e   % -0.4e\n',Mal,Mar);   
%*   fprintf('\n');
  
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1 == 0

  %----------------------------------------------------------------------
  % magnetosonic shocks: DW method
  %----------------------------------------------------------------------  
  
  rfac = 0.24;				% relaxation factor
  
  fprintf('\n');
  fprintf(' Magnetosonic Waves Speeds:\n\n');  
  fprintf('  iter     wmsl     wmsr \n');
  fprintf(' -------------------------\n');
  
  % iterate to find root
  change = 1 + tol;
  kiter = 0;
  while (tol <= change) && (kiter <= maxit)  

    kiter = kiter + 1;
  
    % upstream variable wave traveling to the left/right
    pt_u = ptstar;
  
    % downstream variables wave traveling to the left
    d_d = pstates(1,1);
    pg_d = pstates(5,1);  
    by_d = pstates(6,1);  
    bz_d = pstates(7,1);    

    ms_shock_speed_mhd;
    wmsl = -wms;
    Mal = abs(wmsl)/Cms;
    C01 = C0;
    Cperp1 = Cperp;  
    Cms1 = wmsl; 
  
    % downstream variables wave traveling to the right
    d_d = pstates(1,end);
    pg_d = pstates(5,end);  
    by_d = pstates(6,end);  
    bz_d = pstates(7,end);    
  
    ms_shock_speed_mhd;
    wmsr = wms;
    Mar = abs(wmsr)/Cms;
    C08 = C0;
    Cperp8 = Cperp;  
    Cms8 = wmsr;   
  
    %----------------------------------------------------------------------
    % magnetosonic jump conditions
    %----------------------------------------------------------------------  
  
    % calculate state left of tangential discontinuity (region 2)
    states(:,2) = jump_ms(wmsl,ptstar,ptl,states(:,1));
    pstates(:,2) = cons2prim_mhd_1d(nunks,1,gamma,states(:,2),bx);

    % calculate state right of tangential discontinuity (region 3)
    states(:,3) = jump_ms(wmsr,ptstar,ptr,states(:,4));
    pstates(:,3) = cons2prim_mhd_1d(nunks,1,gamma,states(:,3),bx);

%*     if kiter == 1
%*       fprintf('  %04d   %0.4f   %0.4f\n',kiter,wmsl,wmsr);     
%*     end
%* 
%*     if mod(kiter,100) == 0
      fprintf('  %04d   %0.4f   %0.4f\n',kiter,wmsl,wmsr);       
%*     end


    %----------------------------------------------------------------------
    % Actual iteration step
    %----------------------------------------------------------------------  

    % x-velocity in regions 4 and 5.
    vx4 = pstates(2,2);
    vx5 = pstates(2,3);  
  
    dvx = vx5 - vx4;
  
    DJ4 = ms_jacob_p(gamma,wmsl,states(:,1),states(:,2));
    DJ5 = ms_jacob_p(gamma,wmsr,states(:,4),states(:,3));  
    DJ = DJ4 - DJ5;

    ptnew = ptstar + dvx/(DJ*rfac);
  
    %calulate relative change and update solution
    change = 2*abs((ptnew - ptstar)/(ptnew + ptstar));
    ptstar = ptnew;

  end
  
%*   fprintf('  %04d   %0.4f   %0.4f\n',kiter,wmsl,wmsr);         
  
  fprintf('\n');  
  fprintf(' Mach Numbers:\n\n');
  fprintf('    Mal       Mar\n');
  fprintf(' --------------------\n');
  fprintf('  %0.4f   %0.4f\n',Mal,Mar);   
  fprintf('\n');

end