%******************************************************************************* 
%* Program: guess_mhd.m
%* Description: Computes initial guess of the state vector q.  Called by
%*              exact_mhd_1d.m, and exact nonlinear Riemann sovler for
%*              the 1D Magnetohydrodynamic equations.
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
%* Magnetosonic shock, Bn = 0, guess the common total pressure in regions
%* 4 and 5.  Total pressure (pg + pb) is unchanged across tangential
%* discontinuity. 
%******************************************************************************* 

if strcmp(init_guess,'user')  
 
  if abs(bx) > 0
    bperp(2) = bp2;
    bperp(4) = bp4;
    bperp(7) = bp7;  
    psi(3) = psi3;
  else
    bperp(2) = bp2;
    bperp(7) = bp7;  
  end

elseif strcmp(init_guess,'hlld')
  adaptive = 'off';
  [qhlld,fhlld,vhlld] = hlld(nunks,gamma,bx,pstates(:,1),pstates(:,end));
  
  pstates(:,2) = qhlld(:,2);
  pstates(:,3) = qhlld(:,3);
  pstates(:,4) = qhlld(:,3);  
  pstates(:,5) = qhlld(:,4);
  pstates(:,6) = qhlld(:,4);
  pstates(:,7) = qhlld(:,5);
  
  bperp(2) = sqrt(pstates(6,2)^2 + pstates(7,2)^2);
  bperp(4) = sqrt(pstates(6,4)^2 + pstates(7,4)^2);  
  bperp(7) = sqrt(pstates(6,7)^2 + pstates(7,7)^2);
  psi(3) = atan2(pstates(7,4),pstates(6,4));
  
elseif strcmp(init_guess,'RJ')  
  % initialize average state variables
  qavg = zeros(neigen+1,1);  

  davg = 0.5*(dl + dr);
  vxavg = 0.5*(vxl + vxr);
  vyavg = 0.5*(vyl + vyr);  
  vzavg = 0.5*(vzl + vzr); 
  ptavg = 0.5*(ptl + ptr);  
  byavg = 0.5*(byl + byr);
  bzavg = 0.5*(bzl + bzr);  
  bpavg = sqrt(byavg^2 + bzavg^2);
  
  qavg(1) = davg;
  qavg(2) = vxavg;
  qavg(3) = vyavg;
  qavg(4) = vzavg;
  qavg(5) = ptavg;
  qavg(6) = byavg;
  qavg(7) = bzavg;
  qavg(8) = bx;      

  % calculate eigenvector and eigenmatrices
  [lambda,ematL,ematR] = eigen_mhd(neigen,gamma,qavg,bpavg);
  
  ematR = ematR';
  ematL = ematL';

% calculate initial guess from eq. 3.19 of [2].
  qguess = zeros(nvars,nregio);
  qguess(:,1) = states(:,1);
  qguess(:,end) = states(:,end);

  alpha = ematL*(states(:,end) - states(:,1));
  for i=2:nregio-1
    tmp = 0;
  for kwave=i:nwaves
    tmp = tmp + alpha(kwave)*ematR(:,kwave);  
  end
    qguess(:,i) = qguess(:,end) - tmp;   
  end  

  pstates(:,:) = cons2prim_mhd_1d(nunks,nregio,gamma,qguess(:,:),bx);  
  
  bperp(2) = sqrt(pstates(6,2)^2 + pstates(7,2)^2);
  bperp(4) = sqrt(pstates(6,4)^2 + pstates(7,4)^2);  
  bperp(7) = sqrt(pstates(6,7)^2 + pstates(7,7)^2);
  psi(3) = atan2(pstates(7,3),pstates(6,3));

  wcond = 0;
  cnt = 0;
  while (wcond < 1) && (cnt < 100)

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

    wspeed_check;

    cnt = cnt + 1;
  end

end
