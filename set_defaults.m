%******************************************************************************* 
%* Program: set_defualts.m
%* Description: Sets defualt values for undefined variables.
%* Author: Andrew Kercher 
%-------------------------------------------------------------------------------
%******************************************************************************* 

% parameters
if ~exist('nunks','var')
  nunks = 7;				% number of unknowns
end
if ~exist('nvars','var')
  nvars = 7;				% number of variables
end
if ~exist('neigen','var')
  neigen = 7;				% number of eigenvectors
end
if ~exist('nwaves','var')
  nwaves = 7;				% number of waves
end
if ~exist('nregio','var')
  nregio = 8;				% number of regions/states
end

% see if values are user defined, otherwise set to default values
% number of iterations, tolerance and relaxation factor
if ~exist('state_dimensions','var')
  state_dimensions = 'none';			% dimensional or nondimensional?
end
if ~exist('maxit','var')
  maxit = 20;
end
if ~exist('tol','var')
  tol = 1.0e-6;				% tolerance
end

if ~exist('rfac_k','var')
  rfac_k = 0.5;				% relaxation factor Newton's method
end
if ~exist('rfac_min','var')
  rfac_min = 1.0e-10;			% minimum value of relaxation factor
end
if ~exist('large_n','var')
  large_n = 1e6;			% large number
end
if ~exist('ndim','var')
  ndim = 1;				% number of dimensions for approximation
end
if ~exist('ncell','var')
  ncell = 128;				% number of cells  
end
if ~exist('Lx','var')
  Lx1 = 1;				% domain length
end
if ~exist('x1max','var')
  x1max = 1;				% right bound
end
if ~exist('x1min','var')
  x1min = 0;				% left bound
end
if ~exist('xd0','var')
  xd0 = 0.5;				% initial position of discontinuity
end
if ~exist('delx','var')
  delx = 1.0e-3;			% step-size for FD approx. of Jacobian
end
if ~exist('method','var')
  method = 'RJ';				% RJ or DW
end
if ~exist('derivative_method','var')
  derivative_method = 'RJ';		% RJ or DW
end
if ~exist('Jacobian_method','var')
  Jacobian_method = 'approx';		% nested or approx
end
if ~exist('sampling','var')
  sampling = 'off';		% sample solution, on or off
end
if ~exist('plot_type','var')
  plot_type = 'exact_struc';			% none, exact_struc, exact_cells, approx, exact_approx
end
if ~exist('plot_sym','var')
  plot_sym = '-';				% symbol
end
if ~exist('LW','var')
  LW = 2;				% LineWidth
end
if ~exist('MS','var')
  MS = 6;				% MarkerSize, if point will not
                                        % change unless > 10
end
if ~exist('error_type','var')
  error_type = 'zone_avg';		% zone_avg or max_change
end
if ~exist('init_guess','var')
  init_guess = 'hlld';			% hlld, RJ or user
end
if ~exist('athena_output','var')
  athena_output = 'none';			% none or save
end
if ~exist('bperp_output','var')
  bperp_output = 'none';			% none or save
end

if ~exist('ndim','var')
  ndim = 1;
end
if ~exist('ntout','var')
  ntout = 1;				% number of outputs
end
if ~exist('flux','var')
  flux = 'hlld';
end
if ~exist('order','var')
  order = '3p';
end

if ~exist('approx_solver','var')
  approx_solver = 'none';
end

if ~exist('viscosity','var')
  viscosity = 'disabled';
end
if ~exist('resistivity','var')
  resistivity = 'disabled';
end
if ~exist('conduction','var')
  conduction = 'disabled';
end

if ~exist('lstate','var')
  lstate = 1;
end
if ~exist('rstate','var')
  rstate = 8;
end

if ~exist('athname','var')
  athname = prob;
end
