%******************************************************************************* 
%* Program: states_dim.m
%* Description: Dimensionalizes initial states.
%* Author: Andrew Kercher 
%* References: 
%-------------------------------------------------------------------------------
%******************************************************************************* 
function [wprim,ucons,bperp,wspd,grid_1d,tf,bn,dim] = ...
      states_dim(nunks,ndof,gamma,ds,pgs,x1min,Lx1,wprim,bperp,wspd,grid_1d,tf,bn)
      
pgs = gamma*pgs;
c0s = sqrt(pgs/ds);
C0s = sqrt(ds*pgs);

% compute dimensional states
wprim(1,:) = wprim(1,:)*dl;
wprim(2,:) = wprim(2,:)*c0s;
wprim(3,:) = wprim(3,:)*c0s;
wprim(4,:) = wprim(4,:)*c0s;
wprim(5,:) = wprim(5,:)*pgs;
wprim(6,:) = wprim(6,:)*sqrt(pgs);
wprim(7,:) = wprim(7,:)*sqrt(pgs);

bperp(:) = sqrt(pgs)*bperp(:);

% compute dimensional wave speeds
wspd = C0s*wspd;

% convert parameters
tf = Lx1*tf/c0s;
bn = bn*sqrt(pgs);

% translate from interval [0 1] to original domain
grid_1d.xmin = Lx1*(grid_1d.xmin + x1min);
grid_1d.xmax = Lx1*(grid_1d.xmax + x1min);
grid_1d.xpos = Lx1*(grid_1d.xpos + x1min);
grid_1d.xd = Lx1*(grid_1d.xd + xmin);

ucons = prim2cons_mhd_1d(nunks,ndof,gamma,wprim,bn);

dim = 'dimensional';