%******************************************************************************* 
%* Program: sample_mhd.m
%* Description: Calculates solution for wave.  Sampling preformed in terms of
%*              speed wsi = x/t.  Called by exact_mhd_1d.m, a nonlinear Riemann
%*              solvler for the 1D ideal MHD equations.
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
%* Variables: wsi = speed at point i. Defined wsi = xi/t.
%*        pstates = primitive variables (d,vx,vy,vz,en,by,bz)
%*         states = conservative variables (d,d*vx,d*vy,d*vz,en,by,bz)
%*         pstate = primitive variables at point xi 
%*                   (di,vxi,vyi,vzi,eni,byi,bzi)
%* 
%******************************************************************************* 

if abs(bx) > 0
  % define vx across the contact discontinuity  
  vxd = 0.5*(pstates(2,4) + pstates(2,5));
  
  %position of contact discontinuity
  xpos(4,1) = tf*vxd + xd0;
  xpos(4,2) = tf*vxd + xd0;
  
  %---------------------------------------------------------------  
  % fast/slow shocks/rarefactions and rotational discontinuities
  %---------------------------------------------------------------  
  if wsi <= vxd
    %---------------------------------------------------------------  
    % sampling point is to the left of contact discontinuity
    %---------------------------------------------------------------  
  
    % left fast wave, wfl is negative
    wfasti = states(2,1) - abs(wfl);  
  
    % left rotational discontinuity, wal is negative
    wai = states(2,2) - abs(wal);

    % left slow wave, wsl is negative
    wslowi = states(2,3) - abs(wsl);

    if pstates(1,2)*wsi >= wai
      %------------------------------------------------------------------  
      % sampled pont is downstream (post-shock) rotational discontinuity  
      %------------------------------------------------------------------  

        % position of rotation if it exists
        xpos(2,1) = tf*wai/pstates(1,2) + xd0;
        xpos(2,2) = tf*wai/pstates(1,2) + xd0;

      if pstates(5,end/2) <= pstates(5,end/2-1) && 1 == 1
        %---------------------------------------------------------------  
        % left slow rarefaction wave.
        %---------------------------------------------------------------  

        % compute speed at head of slow wave, Csh
        dh = pstates(1,end/2-1);
        vxh = pstates(2,end/2-1);
        pgh = pstates(5,end/2-1);
        bperph = bperp(end/2-1);

        [C0h2,Cah2,Cperph2,Csh2,Cfh2] = Lagrangian_wave_speeds_mhd(gamma,...
                                        dh,pgh,bperph,bx);
        Csh = sqrt(Csh2);

        % speed at head of left slow rarefaction wave
        wsh = dh*vxh - Csh;

        % compute speed at tail of slow wave, Cst
        dt = pstates(1,end/2);
        vxt = pstates(2,end/2);
        pgt = pstates(5,end/2);
        bperpt = bperp(end/2);

        [C0t2,Cat2,Cperpt2,Cst2,Cft2] = Lagrangian_wave_speeds_mhd(gamma,...
                                          dt,pgt,bperpt,bx);
        Cst = sqrt(Cst2);

        % speed at tail of left slow rarefaction wave
        wst = dt*vxt - Cst;

        if dh*wsi <= wsh
          % sampled point is upstream of slow rarefaction  
          pstate(:,icell) = pstates(:,end/2-1);
        else

          if dt*wsi >= wst
            % sampled point is downstream of slow rarefaction
            pstate(:,icell) = pstates(:,end/2);  
          else
            % sampled point is in rarefaction fan
            lsr_cells = [lsr_cells;icell];
            pstate(:,icell) = 0;
          end
        end

        %position of left slow rarefaction: head --> tail
        xpos(3,1) = tf*wsh/dh + xd0;
        xpos(3,2) = tf*wst/dt + xd0;

      else
        %---------------------------------------------------------------  
        % left slow shock.
        %---------------------------------------------------------------  
        if pstates(1,end/2-1)*wsi <= wslowi
          % sampled point is downstream of slow shock  
          pstate(:,icell) = pstates(:,end/2-1);
        else
          % sampled point is left of contact discontinuity    
          pstate(:,icell) = pstates(:,end/2);
        end  

        %position of left slow shock
        xpos(3,1) = tf*wslowi/pstates(1,end/2-1) + xd0;
        xpos(3,2) = tf*wslowi/pstates(1,end/2-1) + xd0;
      end % left slow structure

    else
      %---------------------------------------------------------------  
      % sampled pont is upstream (preshock) rotational discontinuity  
      %---------------------------------------------------------------  
      % check if gas pressure decreases post-shock
      if pstates(5,2) <= pstates(5,1) && 1 == 1
        %---------------------------------------------------------------  
        % left fast rarefaction wave.
        %---------------------------------------------------------------  

        % compute speed at head of wave, Cfh
        dh = pstates(1,1);
        vxh = pstates(2,1);
        pgh = pstates(5,1);
        bperph = bperp(1);

        [C0h2,Cah2,Cperph2,Csh2,Cfh2] = Lagrangian_wave_speeds_mhd(gamma,...
                                        dh,pgh,bperph,bx);
        Cfh = sqrt(Cfh2);
  
        % speed at head of left fast rarefaction wave
        wfh = dh*vxh - Cfh;

        % compute speed at tail of wave, Cft
        dt = pstates(1,2);
        vxt = pstates(2,2);
        pgt = pstates(5,2);
        bperpt = bperp(2);

        [C0t2,Cat2,Cperpt2,Cst2,Cft2] = Lagrangian_wave_speeds_mhd(gamma,...
                                          dt,pgt,bperpt,bx);
        Cft = sqrt(Cft2);

        % speed at tail of right fast rarefaction wave
        wft = dt*vxt - Cft;
    
        if dh*wsi <= wfh
          % sampled point is initial left state  
          pstate(:,icell) = pstates(:,1);
        else
 
          if dt*wsi >= wft
            % sampled point is downstream of fast rarefaction
            pstate(:,icell) = pstates(:,2);  
          else      
            % sampled point is in rarefaction fan
            lfr_cells = [lfr_cells;icell];
            pstate(:,icell) = 0;
          end
        end

        %position of left fast rarefaction: head --> tail
        xpos(1,1) = tf*wfh/dh + xd0;
        xpos(1,2) = tf*wft/dt + xd0;

      else
        %---------------------------------------------------------------  
        % left fast shock.
        %---------------------------------------------------------------  
        if pstates(1,1)*wsi <= wfasti
          % sampled point is left data state  
          pstate(:,icell) = pstates(:,1);
        else
          % sampled point is upstream of rotational discontinuity    
          pstate(:,icell) = pstates(:,2);
	end  
	
        %position of left fast shock
        xpos(1,1) = tf*wfasti/pstates(1,1) + xd0;
        xpos(1,2) = tf*wfasti/pstates(1,1) + xd0;
	
      end % left fast structure

    end % left structures
  else
    %---------------------------------------------------------------    
    % sampling point is to the right of contact discontinuity
    %---------------------------------------------------------------  
  
    % right fast wave
    wfasti = states(2,8) + wfr;  
  
    % right rotational discontinuity
    wai = states(2,7) + war;

    % right slow wave
    wslowi = states(2,6) + wsr;

    if pstates(1,end-1)*wsi <= wai
      %------------------------------------------------------------------  
      % sampled pont is downstream (post-shock) rotational discontinuity  
      %------------------------------------------------------------------  

      % check for rotation
        %position of right rotational discontinuity
        xpos(end-1,1) = tf*wai/pstates(1,end-1) + xd0;
        xpos(end-1,2) = tf*wai/pstates(1,end-1) + xd0;

      % check if gas pressure decrease behind shock
      if pstates(5,end/2+1) <= pstates(5,end/2+2) && 1 == 1
        %---------------------------------------------------------------  
        % right slow rarefaction wave.
        %---------------------------------------------------------------  

        % compute speed at head of slow wave, Csh
        dh = pstates(1,end/2+2);
        vxh = pstates(2,end/2+2);
        pgh = pstates(5,end/2+2);
        bperph = bperp(end/2+2);

        [C0h2,Cah2,Cperph2,Csh2,Cfh2] = Lagrangian_wave_speeds_mhd(gamma,...
                                        dh,pgh,bperph,bx);
        Csh = sqrt(Csh2);

        % speed at head of slow right rarefaction wave
        wsh = dh*vxh + Csh;

        % compute speed at tail of slow wave, Cst
        dt = pstates(1,end/2+1);
        vxt = pstates(2,end/2+1);
        pgt = pstates(5,end/2+1);
        bperpt = bperp(end/2+1);

        [C0t2,Cat2,Cperpt2,Cst2,Cft2] = Lagrangian_wave_speeds_mhd(gamma,...
                                          dt,pgt,bperpt,bx);
        Cst = sqrt(Cst2);

        % speed at tail of right fast rarefaction wave
        wst = dt*vxt + Cst;
    
        if dh*wsi >= wsh
          % sampled point is upstream of slow rarefaction  
          pstate(:,icell) = pstates(:,end/2+2);
        else

          if dt*wsi <= wst
            % sampled point is downsteam of slow rarefaction
            pstate(:,icell) = pstates(:,end/2+1);  
          else
            % sampled point is right rarefaction fan
            rsr_cells = [rsr_cells;icell];
            pstate(:,icell) = 0;
          end
        end

        %position of right slow rarefaction: tail --> head
        xpos(end-2,1) = tf*wst/dt + xd0;
        xpos(end-2,2) = tf*wsh/dh + xd0;

      else
        %---------------------------------------------------------------  
        % right slow shock.
        %---------------------------------------------------------------  
        if pstates(1,end/2+2)*wsi >= wslowi
          % sampled point is right of slow shock  
          pstate(:,icell) = pstates(:,end/2+2);
        else
          % sampled point is right of contact discontinuity    
          pstate(:,icell) = pstates(:,end/2+1);
        end  

        %position of right slow shock
        xpos(end-2,1) = tf*wslowi/pstates(1,end/2+2) + xd0;
        xpos(end-2,2) = tf*wslowi/pstates(1,end/2+2) + xd0;

      end
    else
      %---------------------------------------------------------------  
      % sampled piont is upstream (preshock) rotational discontinuity  
      %---------------------------------------------------------------  

      % check if gas pressure decrease behind shock
      if pstates(5,end-1) <= pstates(5,end) && 1 == 1
        %---------------------------------------------------------------  
        % right fast rarefaction wave.
        %---------------------------------------------------------------  

        % compute speed at head of wave, Cfh
        dh = pstates(1,end);
        vxh = pstates(2,end);
        pgh = pstates(5,end);
        bperph = bperp(end);

        [C0h2,Cah2,Cperph2,Csh2,Cfh2] = Lagrangian_wave_speeds_mhd(gamma,...
                                        dh,pgh,bperph,bx);
        Cfh = sqrt(Cfh2);
  
        % speed at head of right fast rarefaction wave
        wfh = dh*vxh + Cfh;

        % compute speed at tail of wave, Cft
        dt = pstates(1,end-1);
        vxt = pstates(2,end-1);
        pgt = pstates(5,end-1);
        bperpt = bperp(end-1);

        [C0t2,Cat2,Cperpt2,Cst2,Cft2] = Lagrangian_wave_speeds_mhd(gamma,...
                                          dt,pgt,bperpt,bx);
        Cft = sqrt(Cft2);

        % speed at tail of right fast rarefaction wave
        wft = dt*vxt + Cft;
    
        if dh*wsi >= wfh
          % sampled point is initial right state  
          pstate(:,icell) = pstates(:,end);
        else
 
          if dt*wsi <= wft
            % sampled point is downstream of fast rarefaction
            pstate(:,icell) = pstates(:,end-1);  
          else      
            % sampled point is right rarefaction fan
            rfr_cells = [rfr_cells;icell];
            pstate(:,icell) = 0;
          end
        end

        %position of right fast rarefaction: tail --> head
        xpos(end,1) = tf*wft/dt + xd0;
        xpos(end,2) = tf*wfh/dh + xd0;

      else
        %---------------------------------------------------------------  
        % right fast shock.
        %---------------------------------------------------------------  
        if pstates(1,end)*wsi >= wfasti
          % sampled point is right data state  
          pstate(:,icell) = pstates(:,end);
        else
          % sampled point is upstream of rotational discontinuity    
          pstate(:,icell) = pstates(:,end-1);
	end  
	
        %position of right fast shock
        xpos(end,1) = tf*wfasti/pstates(1,end) + xd0;
        xpos(end,2) = tf*wfasti/pstates(1,end) + xd0;	
  
      end      
    end
  end

else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% define vn across the tangential discontinuity  
vxd = 0.5*(pstates(2,2) + pstates(2,7));
  
%position of tangential discontinuity
xpos(4,1) = tf*vxd + xd0;
xpos(4,2) = tf*vxd + xd0;
  
%---------------------------------------------------------------  
% Magnetosonic shocks
%---------------------------------------------------------------  
if wsi <= vxd

  %---------------------------------------------------------------    
  % sampling point is to the left of tangential discontinuity
  %---------------------------------------------------------------  
    
  % check if gas pressure decrease behind shock
  if pstates(5,2) <= pstates(5,1)

    % left fast rarefaction wave.
    % compute speed at head of wave, Cmsh
    C0h2 = gamma*pstates(1,1)*pstates(5,1);
    Cperph2 = pstates(1,1)*bperp(1)*bperp(1);
    Cmsh = sqrt(C0h2 + Cperph2);

    % speed at head of left magnetosonic rarefaction wave
    wmsh = pstates(1,1)*pstates(2,1) - Cmsh;

    % compute speed at tail of wave, Cmst
    C0t2 = gamma*pstates(1,2)*pstates(5,2);
    Cperpt2 = pstates(1,2)*bperp(2)*bperp(2);
    Cmst = sqrt(C0t2 + Cperpt2);
  
    % speed at tail of left magnetosonic rarefaction wave
    wmst = pstates(1,2)*pstates(2,2) - Cmst;

    if pstates(1,1)*wsi <= wmsh
      % sampled point is initial left state  
      pstate(:,icell) = pstates(:,1);
  
    else

      if pstates(1,2)*wsi > wmst

        % sampled point is left of tangential discontinuity
        pstate(:,icell) = pstates(:,2);  

      else
        
        % sampled point is left rarefaction fan
        pstate(:,icell) = 0;

        % define invariants and constants
        R0 = pstates(5,1)/pstates(1,1)^gamma;
        R1 = bperp(1)/pstates(1,1);
        K1 = sqrt(gamma*R0);
        K2 = (R1/K1)^2;

        dold = pstate(1,icell-1);

        K3 = wsi - pstates(2,1) - 2*(K1/K2)*(1 + K2*pstates(1,1)^(1/3))^(3/2);
        change = 1 + tol;
        kiter = 0;  
  
        while (tol <= change) && (kiter <= maxit)  
  
          kiter = kiter + 1;

          fold = 2*(K1/K2)*(1+K2*dold^(1/3))^(3/2) ...
                 + sqrt(gamma*R0*dold^(gamma - 1) + R1*R1*dold) ...
                 + K3;

          dfold = K1*(1+K2*dold^(1/3))^(1/2)*dold^(-3/2) ...
                  + (gamma*(gamma-1)*R0*dold^(gamma - 2) + R1*R1) ...
                  /(2*sqrt(gamma*R0*dold^(gamma - 1) + R1*R1*dold));

          dnew = dold - fold/dfold;

          % compute error
          change = 2*abs((dnew - dold)/(dnew + dold));
           
          dold = dnew;

        end

        dlr = dnew;
        pglr = R0*dlr^gamma;
        bylr = cos(psi(2))*R1*dlr;
        bzlr = sin(psi(2))*R1*dlr;
 
        c0lr2 = gamma*pglr/dlr;
        cperplr2 = (bylr*bylr + bzlr*bzlr)/dlr;
        cmslr = sqrt(c0lr2 + cperplr2);

        vxlr = wsi + cmslr;

        pstate(1,icell) = dlr;  
        pstate(2,icell) = vxlr;  
        pstate(5,icell) = pglr;  
        pstate(6,icell) = bylr;  
        pstate(7,icell) = bzlr;  

      end  
    end
  
    %position of left fast rarefaction head and tail
    xpos(1,1) = tf*wmsh/pstates(1,1) + xd0;
    xpos(1,2) = tf*wmst/pstates(1,2) + xd0;
  
  else
    %---------------------------------------------------------------  
    % left magnetosonic shock.
    %---------------------------------------------------------------  
  
    % left magnetosonic shock wave, wmsl is negative
    wmsi = pstates(1,1)*pstates(2,1) - abs(wmsl);

    if pstates(1,1)*wsi <= wmsi
      % sampled pont is initial left state  
      pstate(:,icell) = pstates(:,1);
  
    else
      % sampled point is left of tangential discontinuity
      pstate(:,icell) = pstates(:,2);    
    end
  
    %position of right left shock
    xpos(1,1) = tf*wmsi/pstates(1,1) + xd0;
    xpos(1,2) = tf*wmsi/pstates(1,1) + xd0;	
  end

else  
  %---------------------------------------------------------------    
  % sampling point is to the right of tangential discontinuity
  %---------------------------------------------------------------  
  
  % check if gas pressure decrease behind shock
  if pstates(5,end-1) <= pstates(5,end)
    %---------------------------------------------------------------  
    % right fast rarefaction wave.
    %---------------------------------------------------------------  

    % compute speed at head of wave, Cmsh
    C0h2 = gamma*pstates(1,end)*pstates(5,end);
    Cperph2 = pstates(1,end)*bperp(end)*bperp(end);
    Cmsh = sqrt(C0h2 + Cperph2);
  
    % compute speed at tail of wave, Cmst
    C0t2 = gamma*pstates(1,end-1)*pstates(5,end-1);
    Cperpt2 = pstates(1,end-1)*bperp(end-1)*bperp(end-1);
    Cmst = sqrt(C0t2 + Cperpt2);

    % speed at tail of right magnetosonic rarefaction wave
    wmst = pstates(1,end-1)*pstates(2,end-1) + Cmst;
  

    % speed at head of right magnetosonic rarefaction wave
    wmsh = pstates(1,end)*pstates(2,end) + Cmsh;
    
    if pstates(1,end)*wsi >= wmsh
      % sampled point is initial right state  
      pstate(:,icell) = pstates(:,end);
  
    else

      if pstates(1,end-1)*wsi <= wmst

        % sampled point is right of tangential discontinuity
        pstate(:,icell) = pstates(:,end-1);  

      else
        
        % sampled point is right rarefaction fan        
        pstate(:,icell) = 0;

        % define invariants and constants, see [3].
        R0 = pstates(5,end)/pstates(1,end)^gamma;
        R1 = bperp(end)/pstates(1,end);
        K1 = sqrt(gamma*R0);
        K2 = (R1/K1)^2;

        dold = pstate(1,icell-1);

        K3 = pstates(2,end) - wsi - 2*(K1/K2)*(1 + K2*pstates(1,end)^(1/3))^(3/2);
        change = 1 + tol;
        kiter = 0;  
  
        while (tol <= change) && (kiter <= maxit)  
  
          kiter = kiter + 1;

          fold = 2*(K1/K2)*(1+K2*dold^(1/3))^(3/2) ...
                 + sqrt(gamma*R0*dold^(gamma - 1) + R1*R1*dold) ...
                 + K3;

          dfold = K1*(1+K2*dold^(1/3))^(1/2)*dold^(-3/2) ...
                  + (gamma*(gamma-1)*R0*dold^(gamma - 2) + R1*R1) ...
                  /(2*sqrt(gamma*R0*dold^(gamma - 1) + R1*R1*dold));

          dnew = dold - fold/dfold;

          % compute error
          change = 2*abs((dnew - dold)/(dnew + dold));
           
          dold = dnew;

        end

        drr = dnew;
        pgrr = R0*drr^gamma;
        byrr = cos(psi(end-1))*R1*drr;
        bzrr = sin(psi(end-1))*R1*drr;
 
        c0rr2 = gamma*pgrr/drr;
        cperprr2 = (byrr*byrr + bzrr*bzrr)/drr;
        cmsrr = sqrt(c0rr2 + cperprr2);

        vxrr = wsi - cmsrr;        

        pstate(1,icell) = drr;  
        pstate(2,icell) = vxrr;  
        pstate(5,icell) = pgrr;  
        pstate(6,icell) = byrr;  
        pstate(7,icell) = bzrr;  

      end
  
    end
  
    %position of right fast rarefaction head and tail
    xpos(end,1) = tf*wmst/pstates(1,end-1) + xd0;
    xpos(end,2) = tf*wmsh/pstates(1,end) + xd0;

  else
    %---------------------------------------------------------------  
    % right magnetosonic shock.
    %---------------------------------------------------------------  
  
    wmsi = pstates(1,end)*pstates(2,end) + wmsr;
    
    if pstates(1,end)*wsi >= wmsi
      % sampled ponit is initial right state  
      pstate(:,icell) = pstates(:,end);
  
    else
      % sampled ponit is right of tangential discontinuity
      pstate(:,icell) = pstates(:,end-1);  
    end
  
    %position of right magnetosonic shock
    xpos(end,1) = tf*wmsi/pstates(1,end) + xd0;
    xpos(end,2) = tf*wmsi/pstates(1,end) + xd0;	
  end  
end
end