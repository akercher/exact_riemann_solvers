%******************************************************************************* 
%* Program: plot_mhd.m
%* Description: Plotting routine for solution to exact nonlinear Riemann
%*              sovler for 1D Magnetohydrodynamic equations.
%* Author: Andrew Kercher 
%-------------------------------------------------------------------------------
%******************************************************************************* 

% compute the ylimits
get_ylims;

if strcmp(plot_sym,'.')
  sym_size = 'MarkerSize';
  SS = MS;
else
  sym_size = 'LineWidth';
  SS = LW;  
end

%-------------------------------------------------------------------------------
%* Plot exact solution 
%-------------------------------------------------------------------------------
if strcmp(plot_type,'exact_struc')
  if ~exist('fig_list','var')
    fig_list = 1:9;
  end
  nfig = length(fig_list);
  for iplot = 1:nfig
    ifig = fig_list(iplot);
    plot_profile;
  end
  return;
end

%-------------------------------------------------------------------------------
%* Plot exact solution only across multiple cells
%-------------------------------------------------------------------------------
if strcmp(plot_type,'exact_cells') 

figure(1)
% density
subplot(3,1,1)
plot(xcell,pstate(1,:),plot_sym,sym_size,SS)
ylim(ylim_d)
ylabel('rho')
title(['Problem: ',prob])

% gas pressure
subplot(3,1,2)
plot(xcell,pstate(5,:),plot_sym,sym_size,SS)
ylim(ylim_pg)    
ylabel('p_g');

% energy density
subplot(3,1,3)
plot(xcell,state(5,:),plot_sym,sym_size,SS)
ylim(ylim_en)      
ylabel('E')

figure(2)
% x-velocity
subplot(3,1,1)
plot(xcell,pstate(2,:),plot_sym,sym_size,SS)
ylim(ylim_vx)  
ylabel('v_x')
title(['Problem: ',prob])

% y-velocity
subplot(3,1,2)
plot(xcell,pstate(3,:),plot_sym,sym_size,SS)
ylim(ylim_vy)
ylabel('v_y')

% z-velocity
subplot(3,1,3)
plot(xcell,pstate(4,:),plot_sym,sym_size,SS)
ylim(ylim_vz)      
ylabel('v_z')

figure(3)
% y-magnetic field
subplot(3,1,1)
plot(xcell,pstate(6,:),plot_sym,sym_size,SS)
ylim(ylim_by)  
ylabel('B_y')
title(['Problem: ',prob])

% z-magnetic field
subplot(3,1,2)
plot(xcell,pstate(7,:),plot_sym,sym_size,SS)
ylim(ylim_bz)
ylabel('B_z')

% rotation angle
subplot(3,1,3)
plot(xcell,psi_x,plot_sym,sym_size,SS)
ylim(ylim_psi)
ylabel('psi')

end

%-------------------------------------------------------------------------------
%* Plot exact solution only at position of structures
%-------------------------------------------------------------------------------
if strcmp(plot_type,'exact_pos') 

  xsol = zeros(2*neigen,1);  
  psol = zeros(nunks,2*nregio);
  esol = zeros(2*nregio,1);  

  for i=1:neigen
    xsol(2*i-1) = xpos(i,1);
    xsol(2*i) = xpos(i,2);  
  end
  
  xsol = [0;xsol;1];
  
  for i=1:nregio
    psol(:,2*i-1) = pstates(:,i);
    psol(:,2*i) = pstates(:,i);  
    esol(2*i-1) = states(5,i);    
    esol(2*i) = states(5,i);      
  end
  
figure(1)
% density
subplot(3,1,1)
plot(xsol,psol(1,:),plot_sym,sym_size,SS)
ylim(ylim_d)
ylabel('rho')
title(['Problem: ',prob])

% gas pressure
subplot(3,1,2)
plot(xsol,psol(5,:),plot_sym,sym_size,SS)
ylim(ylim_pg)    
ylabel('p_g');

% energy density
subplot(3,1,3)
plot(xsol,esol,plot_sym,sym_size,SS)
ylim(ylim_en)      
ylabel('en')

figure(2)
% x-velocity
subplot(3,1,1)
plot(xsol,psol(2,:),plot_sym,sym_size,SS)
ylim(ylim_vx)  
ylabel('v_x')
title(['Problem: ',prob])

% y-velocity
subplot(3,1,2)
plot(xsol,psol(3,:),plot_sym,sym_size,SS)
ylim(ylim_vy)
ylabel('v_y')

% z-velocity
subplot(3,1,3)
plot(xsol,psol(4,:),plot_sym,sym_size,SS)
ylim(ylim_vz)      
ylabel('v_z')

figure(3)
% y-magnetic field
subplot(3,1,1)
plot(xsol,psol(6,:),plot_sym,sym_size,SS)
ylim(ylim_by)  
ylabel('B_y')
title(['Problem: ',prob])

% z-magnetic field
subplot(3,1,2)
plot(xsol,psol(7,:),plot_sym,sym_size,SS)
ylim(ylim_bz)
ylabel('B_z')

% rotation angle
subplot(3,1,3)
plot(xsol,(atan2(psol(7,:),psol(6,:))*180/pi),plot_sym,sym_size,SS)
ylim(ylim_psi)
ylabel('psi')

end

%-------------------------------------------------------------------------------
%* Plot approximate solutions
%-------------------------------------------------------------------------------
if strcmp(plot_type,'approx')

  LW = 16;
  FS = 32;
  FSxy = 24;
  
  % include directory with approximate Riemann solver.
  addpath('./src/1d_mhd_approx');
  
  % approximate Riemann solver
%*   mhd1Dcc;

xtick_labels = [{'0.5';'0.6';'0.7'}];
xtick_val = [0.5;0.6;0.7];		  

ytick_labels = [{'0.3';'0.6';}];
ytick_val = [0.3;0.6];		  

% density
figure(1)
clf
plot(xcell,pstate(1,:),'ob','LineWidth',LW)
hold on
plot(xcell,rho(niter,:),'-r','LineWidth',LW)
xlim([0.45 0.7])
ylim([0.05 0.85])
xlabel('x','FontSize',FS);
ylabel('\rho','FontSize',FS);
h_leg = legend('exact','approx',1);
set(h_leg,'FontSize',FS);

set(gca,'XTick',xtick_val);
set(gca,'XTickLabel',xtick_labels);
set(gca,'YTick',ytick_val);
set(gca,'YTickLabel',ytick_labels);
set(gca,'FontSize',FSxy);

com = sprintf('print -depsc %s_approx.eps',prob);
eval(com)
  
%*   figure
%*   % density
%*   subplot(3,1,1)
%*   plot(xcell,pstate(1,:),'o')
%*   hold on
%*   plot(xc,rho(niter,:),'xr')  
%*   ylim(ylim_d) 
%*   ylabel('rho')  
%*   title(['Density'])
%*   legend('Exact','Aprrox.',3)

end

%-------------------------------------------------------------------------------
%* Plot exact and approximate solutions
%-------------------------------------------------------------------------------
if strcmp(plot_type,'exact_approx')
  
  % include directory with approximate Riemann solver.
  addpath('./src/1d_mhd_approx');

  % approximate Riemann solver
  mhd1Dcc;
  
  figure
  % density
  subplot(3,1,1)
  plot(xcell,pstate(1,:),'o')
  hold on
  plot(xc,rho(niter,:),'xr')  
  ylim(ylim_d) 
  ylabel('rho')  
%*   title(['Density'])
%*   legend('Exact','Aprrox.',3)

  % gas pressure
  subplot(3,1,2)
  plot(xcell,pstate(5,:),'o')
  hold on
  plot(xc,p(niter,:),'xr')
  ylim(ylim_pg)  
  ylabel('p_g')
%*   title(['Pressure'])
%*   legend('Exact','Aprrox.',2)
  
  % energy density
  subplot(3,1,3)
  plot(xcell,state(5,:),'o')
  hold on
  plot(xc,en(niter,:),'xr')
  ylim(ylim_en)
  ylabel('en')
%*   title(['Energy Density'])
%*   legend('Exact','Aprrox.',3)  
  
%*   if ncell < 30
%*     print -dpng magnetsonic_1_low.png  
%*   else
%*     print -dpng magnetsonic_1_high.png    
%*   end
  
  figure
  % x-velocity
  subplot(3,1,1)
  plot(xcell,pstate(2,:),'o')
  hold on
  plot(xc,vx(niter,:),'xr')
  ylim(ylim_vx)  
  ylabel('v_x')
%*   title(['x-velocity'])
%*   legend('Exact','Aprrox.',3)  
  
  % y-velocity
  subplot(3,1,2)
  plot(xcell,pstate(3,:),'o')
  hold on
  plot(xc,vy(niter,:),'xr')
  ylim(ylim_vy)
  ylabel('v_y')
%*   title(['y-velocity'])
%*   legend('Exact','Aprrox.',3)  

  % z-velocity
  subplot(3,1,3)
  plot(xcell,pstate(4,:),'o')
  hold on
  plot(xc,vz(niter,:),'xr')
  ylim(ylim_vz)  
  ylabel('v_z')
%*   title(['z-velocity'])
%*   legend('Exact','Aprrox.',3)  

  
%*   if ncell < 30
%*     print -dpng magnetsonic_2_low.png  
%*   else
%*     print -dpng magnetsonic_2_high.png    
%*   end
  
  figure
  % y-magnetic field  
  subplot(3,1,1)
  plot(xcell,pstate(6,:),'o')
  hold on  
  plot(xc,By(niter,:),'xr')
  ylim(ylim_by)  
  ylabel('B_y')
%*   title(['y-magnetic field'])
%*   legend('Exact','Aprrox.',2)  

  % z-magnetic field  
  subplot(3,1,2)
  plot(xcell,pstate(7,:),'o')
  hold on  
  plot(xc,Bz(niter,:),'xr')
  ylim(ylim_bz)  
  ylabel('B_z')
%*   title(['z-magnetic field'])
%*   legend('Exact','Aprrox.',2)  

  % rotation angle in degrees 
  subplot(3,1,3)
  plot(xcell,psi_x,'o')
  hold on  
  plot(xc,180*(atan(Bz(niter,:)./By(niter,:)))/pi,'xr')
  ylim(ylim_psi)  
  ylabel('psi')
%*   title(['Rotation Angle'])
%*   legend('Exact','Aprrox.',3)  

%*   if ncell < 30
%*     print -dpng magnetsonic_3_low.png  
%*   else
%*     print -dpng magnetsonic_3_high.png    
%*   end
  
end
