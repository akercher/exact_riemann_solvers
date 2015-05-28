%******************************************************************************* 
%* Program: states_nondim.m
%* Description: Nondimensionalizes initial states.
%* Author: Andrew Kercher 
%* References: 
%-------------------------------------------------------------------------------
%******************************************************************************* 
function [wprim,ucons,grid_1d,tf,bn,dim] = states_nondim(nunks,ndof,gamma,wprim,grid_1d,tf,bn)

% define scaling factors
istate = 1;  
if wprim(1,end) > wprim(1,1)
  istate = 8;
end  
ds = wprim(1,istate);
pgs = gamma*wprim(5,istate);
c0s = sqrt(pgs/ds);
xmin = grid_1d.xmin;
xmax = grid_1d.xmax;
xd = grid_1d.xd;
Lx = abs(xmax - xmin);

% define non-dimensional states
wprim(1,:) = wprim(1,:)/ds;
wprim(2,:) = wprim(2,:)/c0s;
wprim(3,:) = wprim(3,:)/c0s;
wprim(4,:) = wprim(4,:)/c0s;
wprim(5,:) = wprim(5,:)/(pgs);
wprim(6,:) = wprim(6,:)/sqrt(pgs);
wprim(7,:) = wprim(7,:)/sqrt(pgs);

% convert parameters
tf = c0s*tf/Lx;
bn = bn/sqrt(pgs);

% translate to interval [0 1]
grid_1d.xmin = (grid_1d.xmin - xmin)/Lx;
grid_1d.xmax = (grid_1d.xmax - xmin)/Lx;
grid_1d.xd = (grid_1d.xd - xmin)/Lx;

ucons = prim2cons_mhd_1d(nunks,ndof,gamma,wprim,bn);

dim = 'none';