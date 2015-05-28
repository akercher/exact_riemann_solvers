%******************************************************************************* 
%* Program: exact_mhd_run.m
%* Description: Start up file to run a nonlinear Riemann sovler for
%*              1D Magnetohydrodynamic equations.
%* Author: Andrew Kercher 
%-------------------------------------------------------------------------------
%******************************************************************************* 

if ~exist('run_type','var')
   run_type = 'new';			% new or continued    
end  
 
if strcmp(run_type,'new')
  var_list = {prob};
  clear_all_except;
  prob = var_list{1};

  load_paths;
  pwd_top = sprintf('%s',pwd);

  run_type = 'new';			% new or continued  

  % run configuration file
  init_prob = sprintf('init_%s',prob);
  eval(init_prob);

  % set default values

  set_defaults;
  
if (lstate > 1) || (rstate < 8)
    prob_data = sprintf('%s/dat/%s.dat',pwd,prob);
    load(prob_data,'pstates','states','bperp','psi','bx','grid_1d','state_dimensions');  

    % nondimensionalize final time using original variables
    if ~strcmp(state_dimensions,'none')
      Lx = abs(x1max - x1min);
      c0s = sqrt(gamma*pgl/dl);
      if dr > dl
        c0s = sqrt(gamma*pgr/dr);  
      end
      tf = c0s*tf/Lx;
    end
  
%*     state_dimensions = 'dimensional';
    dl = pstates(1,lstate);
    vxl = pstates(2,lstate);
    vyl = pstates(3,lstate);
    vzl = pstates(4,lstate);
    pgl = pstates(5,lstate);
    byl = pstates(6,lstate);
    bzl = pstates(7,lstate);
    kel = 0.5*(vxl*vxl + vyl*vyl + vzl*vzl);
    pbl = 0.5*(bx*bx + byl*byl + bzl*bzl);
    enl = pgl/(gamma - 1) + dl*kel + pbl;
    ptl = pgl + pbl;

    % initial right state  
    dr = pstates(1,rstate);
    vxr = pstates(2,rstate);
    vyr = pstates(3,rstate);
    vzr = pstates(4,rstate);
    pgr = pstates(5,rstate);
    byr = pstates(6,rstate);
    bzr = pstates(7,rstate);
    ker = 0.5*(vxr*vxr + vyr*vyr + vzr*vzr);
    pbr = 0.5*(bx*bx + byr*byr + bzr*bzr);
    enr = pgr/(gamma - 1) + dr*ker + pbr;
    ptr = pgr + pbr;
  
    if mod(lstate,2) == 0
      bp2 = bperp(lstate);
      if lstate > 2
        bp4 = bperp(lstate);
        psi3 = psi(lstate);
      else
        psi3 = psi(lstate+1);
        if rstate > 4
          bp4 = bperp(4);
        else 
          bp4 = bperp(rstate);
        end
      end
    else
      bp2 = bperp(lstate+1);  
      if lstate > 1
        bp4 = bperp(lstate+1);    
        psi3 = psi(lstate);
      else
        if rstate >= 4
          bp4 = bperp(4);
          psi3 = psi(3);  
        else 
          bp4 = bperp(rstate);
          psi3 = psi(rstate);  
        end      
      end
    end
  
    if rstate > 7
      bp7 = bperp(rstate-1);
    else
      bp7 = bperp(rstate);  
    end

    init_guess = 'user';

%* keyboard
    end
    
else
  % set default values
  set_defaults;
  
  dl = pstates(1,lstate);
  vxl = pstates(2,lstate);
  vyl = pstates(3,lstate);
  vzl = pstates(4,lstate);
  pgl = pstates(5,lstate);
  byl = pstates(6,lstate);
  bzl = pstates(7,lstate);
  kel = 0.5*(vxl*vxl + vyl*vyl + vzl*vzl);
  pbl = 0.5*(bx*bx + byl*byl + bzl*bzl);
  enl = pgl/(gamma - 1) + dl*kel + pbl;
  ptl = pgl + pbl;

  % initial right state  
  dr = pstates(1,rstate);
  vxr = pstates(2,rstate);
  vyr = pstates(3,rstate);
  vzr = pstates(4,rstate);
  pgr = pstates(5,rstate);
  byr = pstates(6,rstate);
  bzr = pstates(7,rstate);
  ker = 0.5*(vxr*vxr + vyr*vyr + vzr*vzr);
  pbr = 0.5*(bx*bx + byr*byr + bzr*bzr);
  enr = pgr/(gamma - 1) + dr*ker + pbr;
  ptr = pgr + pbr;
    
  init_guess = 'user';			% user or RJ  
  
end

% precision for long format
fmt_pre = abs(log10(tol));

% execute main file
exact_mhd_main;  

for j=1:nregio
  if abs(psi(j)) > pi
    fprintf(' WARNING: psi in region %d > pi\n',j);
  end
end

prob_data = sprintf('%s/dat/%s.dat',pwd,prob);
if (lstate == 1) && (rstate == 8)
  save(prob_data);
end
  
