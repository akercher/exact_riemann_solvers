%******************************************************************************* 
%* Program: ms_jacob_fd.m
%* Description: Approximates Jacobian with finite differences.
%* Author: Andrew Kercher 
%* References: 
%-------------------------------------------------------------------------------
%******************************************************************************* 

nx = length(xvec);

DJ = zeros(nx,nx);
fxx = zeros(nx,1);

% save current state variables
pstates_old = pstates;

for j=1:nx
  
  xx = xvec;
  
  xx(j) = xvec(j) + delx;
  
  bperp(2:4) = xx(1);
  bperp(5:7) = xx(2);
  
  % compute states based on updated independent variables
  ms_rstates_mhd;

  % difference in state variables, we want this to converge to zero
  fxx(1) = pstates(2,end-1) - pstates(2,2);
  fxx(2) = pstates(5,end-1) - pstates(5,2) + 0.5*(bperp(end-1)^2 -  bperp(2)^2);
  
  DJ(:,j) = (fxx - ftd)/delx;

end  

% revert to previous state variables
pstates = pstates_old;
bperp(2:4) = xvec(1);
bperp(5:7) = xvec(2);
