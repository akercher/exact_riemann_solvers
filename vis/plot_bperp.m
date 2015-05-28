%******************************************************************************* 
%* Program: plot_bperp.m
%* Description: Plotting routine for perpendicular magnetic field
%* Author: Andrew Kercher 
%-------------------------------------------------------------------------------
%******************************************************************************* 

%-------------------------------------------------------------------------------
%* Plots evolution of perpendicular magnetic field
%-------------------------------------------------------------------------------

npsi = 32;
SS = 16;

figure(4)
clf
hold on
% plot transitions
for iregio=2:nregio
  psi0 = min(psi(iregio-1),psi(iregio));
  psi1 = max(psi(iregio-1),psi(iregio));
  bperp0 = bperp(iregio-1);  
  bperp1 = bperp(iregio);    
  
  dpsi = abs(psi1 - psi0);
  if dpsi < tol
    h1 = plot(pstates(6,iregio-1:iregio),pstates(7,iregio-1:iregio),'--r');
    set(h1,'LineWidth',LW);  
  else  
    theta = psi0:dpsi/(npsi-1):psi1;
    [by,bz] = pol2cart(theta,bperp(iregio));
    h1 = plot(by,bz,'--r');    
    set(h1,'LineWidth',LW);

  end
end  

% plot states
for iregio=1:nregio
  h1 = plot(pstates(6,iregio),pstates(7,iregio),'k.');
  set(h1,'MarkerSize',SS);
end

% set axis
bymin = min(pstates(6,:));
bymax = max(pstates(6,:));
bzmin = min(pstates(7,:));
bzmax = max(pstates(7,:));

xlim1 = 0;
xlim2 = 1;
ylim1 = 0;
ylim2 = 1;
if bymin < 0
  xlim1 = -1;
end
if bymax < 0
  xlim2 = 0;
end
if bzmin < 0
  ylim1 = -1;
end
if bzmax < 0
  ylim2 = 0;
end

xlim([xlim1 xlim2])
ylim([ylim1 ylim2])

title(['problem= ',prob]);

if strcmp(bperp_output,'save')
  figid = sprintf('%bperp_s_%s_%s_%d',prob,flux,order,i);
  com = sprintf('print -depsc %s/fig/%s.eps',pwd_top,figid);
  eval(com);
end
