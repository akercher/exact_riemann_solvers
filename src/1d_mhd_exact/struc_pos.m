%******************************************************************************* 
%* Program: struc_pos.m
%* Description: Calculates the postions of the structure (wave/shock).
%* Author: Andrew Kercher 
%* References: 
%*         [1] Toro, E. F., "Riemann Solvers and Numerical Methods for
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

function [xpos,vspd] = struc_pos(gamma,pstates,bperp,wspd,xd0,tf,bn);

xpos = zeros(length(wspd)+1,1);  
vspd = zeros(length(wspd)+1,1);  
wfl = wspd(1);
wal = wspd(2);
wsl = wspd(3);
wsr = wspd(4);
war = wspd(5);
wfr = wspd(6);

if abs(bn) > 0

  % define normal velocity at contact discontinuity  
  vcd = 0.5*(pstates(2,4) + pstates(2,5));
  
  %position of contact discontinuity
  xpos(4,1) = tf*vcd + xd0;
  xpos(4,2) = tf*vcd + xd0;

  % velocities for left (minus) fast, Alfven and slow.
  % Left wave speeds are negative.
  vfl = pstates(2,1) - abs(wfl)/pstates(1,1);
  val = pstates(2,2) - abs(wal)/pstates(1,2);
  vsl = pstates(2,3) - abs(wsl)/pstates(1,3);  

  % position of fast slow shock if it exists  
  xpos(1,:) = tf*vfl + xd0;
  
  % position of left Alfven wave if it exists
  xpos(2,:) = tf*val + xd0;

  % position of left slow shock if it exists  
  xpos(3,:) = tf*vsl + xd0;
  
  % velocities for right (plus) fast, Alfven and slow.
  vfr = pstates(2,end) + wfr/pstates(1,end);
  var = pstates(2,end-1) + war/pstates(1,end-1);
  vsr = pstates(2,end-2) + wsr/pstates(1,end-2);

  % position of right fast shock if it exists  
  xpos(end,:) = tf*vfr + xd0;
  
  % position of right Alfven wave if it exists
  xpos(end-1,:) = tf*var + xd0;

  % position of right slow shock if it exists  
  xpos(end-2,:) = tf*vsr + xd0;
  
  vspd(1) = vfl;
  vspd(2) = val;  
  vspd(3) = vsl;    
  vspd(4) = vcd;      
  vspd(5) = vsr;      
  vspd(6) = var;      
  vspd(7) = vfr;        
  
  if pstates(5,end/2) <= pstates(5,end/2-1)
    %---------------------------------------------------------------  
    % left slow rarefaction wave.
    %---------------------------------------------------------------  
  
    % compute speed at head of slow wave, Csh
    [C0h2,Cah2,Cperph2,Csh2,Cfh2] = Lagrangian_wave_speeds_mhd(gamma,...
                                     pstates(1,end/2-1),pstates(5,end/2-1),bperp(end/2-1),bn);

    % speed at head of left slow rarefaction wave
    vsh = pstates(2,end/2-1) - sqrt(Csh2)/pstates(1,end/2-1);

    % compute speed at tail of slow wave, Cst
    [C0t2,Cat2,Cperpt2,Cst2,Cft2] = Lagrangian_wave_speeds_mhd(gamma,...
                                      pstates(1,end/2),pstates(5,end/2),bperp(end/2),bn);

    % speed at tail of left slow rarefaction wave
    vst = pstates(2,end/2) - sqrt(Cst2)/pstates(1,end/2);

    % position of left slow rarefaction head --> tail
    xpos(3,1) = tf*vsh + xd0;
    xpos(3,2) = tf*vst + xd0;
  
  end % end left slow rarefaction
  
  if pstates(5,2) <= pstates(5,1)
    %---------------------------------------------------------------  
    % left fast rarefaction wave.
    %---------------------------------------------------------------  

    % compute speed at head of wave, Cfh
    [C0h2,Cah2,Cperph2,Csh2,Cfh2] = Lagrangian_wave_speeds_mhd(gamma,...
                                      pstates(1,1),pstates(5,1),bperp(1),bn);
    Cfh = sqrt(Cfh2);
  
    % speed at head of right fast rarefaction wave
    vfh = pstates(2,1) - sqrt(Cfh2)/pstates(1,1);
  
    % compute speed at tail of wave, Cft
    [C0t2,Cat2,Cperpt2,Cst2,Cft2] = Lagrangian_wave_speeds_mhd(gamma,...
                                      pstates(1,2),pstates(5,2),bperp(2),bn);
  
    % speed at tail of left fast rarefaction wave
    vft = pstates(2,2) - sqrt(Cft2)/pstates(1,2);
  
    % position of left fast rarefaction head --> tail
    xpos(1,1) = tf*vfh + xd0;
    xpos(1,2) = tf*vft + xd0;

  end % left fast structure
  
  if pstates(5,end-3) <= pstates(5,end-2)
    %---------------------------------------------------------------  
    % right slow rarefaction wave.
    %---------------------------------------------------------------  
  
    % compute speed at head of slow wave, Csh
    [C0h2,Cah2,Cperph2,Csh2,Cfh2] = Lagrangian_wave_speeds_mhd(gamma,...
                                      pstates(1,end-2),pstates(5,end-2),bperp(end-2),bn);

    % speed at head of right slow rarefaction wave
    vsh = pstates(2,end-2) + sqrt(Csh2)/pstates(1,end-2);

    % compute speed at tail of slow wave, Cst
    [C0t2,Cat2,Cperpt2,Cst2,Cft2] = Lagrangian_wave_speeds_mhd(gamma,...
                                      pstates(1,end-3),pstates(5,end-3),bperp(end-3),bn);

    % speed at tail of right slow rarefaction wave
    vst = pstates(2,end-3) + sqrt(Cst2)/pstates(1,end-3);

    % position of right slow rarefaction: tail --> head
    xpos(end-2,1) = tf*vst + xd0;
    xpos(end-2,2) = tf*vsh + xd0;
  
  end % end right slow rarefaction
  
  if pstates(5,end-1) <= pstates(5,end)
    %---------------------------------------------------------------  
    % right fast rarefaction wave.
    %---------------------------------------------------------------  

    % compute speed at head of wave, Cfh
    [C0h2,Cah2,Cperph2,Csh2,Cfh2] = Lagrangian_wave_speeds_mhd(gamma,...
                                      pstates(1,end),pstates(5,end),bperp(end),bn);

    % speed at head of right fast rarefaction wave
    vfh = pstates(2,end) + sqrt(Cfh2)/pstates(1,end);
  
    % compute speed at tail of wave, Cft
    [C0t2,Cat2,Cperpt2,Cst2,Cft2] = Lagrangian_wave_speeds_mhd(gamma,...
                                      pstates(1,end-1),pstates(5,end-1),bperp(end-1),bn);
  
    % speed at tail of right fast rarefaction wave
    vft = pstates(2,end-1) + sqrt(Cft2)/pstates(1,end-1);
  
    % position of left fast rarefaction tail --> head
    xpos(end,1) = tf*vft + xd0;
    xpos(end,2) = tf*vfh + xd0;

  end % right fast structure
    
  
  
end  