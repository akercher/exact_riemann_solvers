%******************************************************************************* 
%* Program: get_ylims.m
%* Description: Computes ylimit for plotting. 
%* Author: Andrew Kercher 
%* References: 
%-------------------------------------------------------------------------------
%******************************************************************************* 

if ~exist('xstart','var')
  xstart = 0;				% starting position for plotting
end

if ~exist('xend','var')
  xend = 1;				% ending position for plotting
end

% calc. for plotting limits
dmin = min(pstates(1,ceil(xstart*end) + 1:floor(xend*end)));
dmax = max(pstates(1,ceil(xstart*end) + 1:floor(xend*end)));
dd = abs(dmax - dmin);
if dd < tol 
  dd = 1;
end  
vxmin = min(pstates(2,ceil(xstart*end) + 1:floor(xend*end)));
vxmax = max(pstates(2,ceil(xstart*end) + 1:floor(xend*end)));
dvx = abs(vxmax - vxmin);
vymin = min(pstates(3,ceil(xstart*end) + 1:floor(xend*end)));
vymax = max(pstates(3,ceil(xstart*end) + 1:floor(xend*end)));
dvy = abs(vymax - vymin);
vzmin = min(pstates(4,ceil(xstart*end) + 1:floor(xend*end)));
vzmax = max(pstates(4,ceil(xstart*end) + 1:floor(xend*end)));
dvz = abs(vzmax - vzmin);
pgmin = min(pstates(5,ceil(xstart*end) + 1:floor(xend*end)));
pgmax = max(pstates(5,ceil(xstart*end) + 1:floor(xend*end)));
dpg = abs(pgmax - pgmin);
enmin = min(states(5,ceil(xstart*end) + 1:floor(xend*end)));
enmax = max(states(5,ceil(xstart*end) + 1:floor(xend*end)));
den = abs(enmax - enmin);
bymin = min(pstates(6,ceil(xstart*end) + 1:floor(xend*end)));
bymax = max(pstates(6,ceil(xstart*end) + 1:floor(xend*end)));
dby = abs(bymax - bymin);
bzmin = min(pstates(7,ceil(xstart*end) + 1:floor(xend*end)));
bzmax = max(pstates(7,ceil(xstart*end) + 1:floor(xend*end)));
dbz = abs(bzmax - bzmin);
psimin = min(psi(ceil(xstart*end) + 1:floor(xend*end)));
psimax = max(psi(ceil(xstart*end) + 1:floor(xend*end)));
dpsi = abs(psimax - psimin);

% y-axis limits for plotting
ylim_d = [dmin - 0.1*dd - 0.1 dmax + 0.1*dd];  
ylim_vx = [vxmin - 0.1*dvx - 0.1 vxmax + 0.1*dvx + 0.1];  
ylim_vy = [vymin - 0.1*dvy - 0.1 vymax + 0.1*dvy + 0.1];  
ylim_vz = [vzmin - 0.1*dvz - 0.1 vzmax + 0.1*dvz + 0.1];  
ylim_pg = [pgmin - 0.1*dpg - 0.1 pgmax + 0.1*dpg + 0.1];  
ylim_en = [enmin - 0.1*den - 0.1 enmax + 0.1*den + 0.1];  
ylim_by = [bymin - 0.1*dby - 0.1 bymax + 0.1*dby + 0.1];  
ylim_bz = [bzmin - 0.1*dbz - 0.1 bzmax + 0.1*dbz + 0.1];  
ylim_psi = [psimin - 0.1*dpsi - 0.1 psimax + 0.1*dpsi + 0.1];  

limits_y = [ylim_d; ylim_vx; ylim_vy; ylim_vz; ylim_pg; ylim_by; ylim_bz; ...
	  ylim_en; ylim_psi];
    
labels_y = {'rho'; 'v_x'; 'v_y'; 'v_z'; 'p_g'; 'b_y'; 'b_z'; 'E'; 'psi'};
if strcmp(athena_output,'latex')
  labels_y = {'\LARGE $\rho$'; '\Large $v_x$'; '\Large $v_y$'; '\Large $v_z$'; '\Large $p_g$'; '\LARGE $B_y$'; '\Large $B_z$'; '\Large $E$'; '$\Huge \psi$'};  
end  
