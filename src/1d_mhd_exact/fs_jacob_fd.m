%******************************************************************************* 
%* Program: fs_jacob_fd.m
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
  
  bperp(2:3) = xx(1);
  bperp(4:5) = xx(2);
  bperp(6:7) = xx(3);
  psi(3:6) = xx(4);

  % check for switch-structures
  if abs(bperp(1)) < tol
    psi(1:2) = psi(3);
  end
  
  if abs(bperp(end)) < tol  
    psi(end-1:end) = psi(3);  
  end
  
  % compute states based on updated independent variables
  fs_rstates_mhd;

  % difference in state variables, we want this to converge to zero
  fxx(1) = pstates(2,5) - pstates(2,4);
  fxx(2) = pstates(3,5) - pstates(3,4);
  fxx(3) = pstates(4,5) - pstates(4,4);
  fxx(4) = pstates(5,5) - pstates(5,4) + 0.5*(bperp(5)^2 -  bperp(4)^2);
  
  DJ(:,j) = (fxx - fcd)/delx;

end  

%* DJ
%* xfdfga


xvec(1:3) = abs(xvec(1:3));
%* if abs(xvec(4) > pi)
%*   xvec(4) = 2*pi - abs(xvec(4));
%* end

% revert to previous state variables
pstates = pstates_old;
bperp(2:3) = xvec(1);
bperp(4:5) = xvec(2);
bperp(6:7) = xvec(3);
psi(3:6) = xvec(4);
