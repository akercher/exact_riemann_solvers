%******************************************************************************* 
%* Program: plot_profile.m
%* Description: plots solution profile.
%* Author: Andrew Kercher 
%* References: 
%-------------------------------------------------------------------------------
%******************************************************************************* 

if strcmp(athena_output,'latex')
grid_1d.xmin = xstart;
grid_1d.xmax = xend;
end

if ifig < 8
  figure(ifig)
%*   clf
    if (xr(1,2) >= grid_1d.xmin) && (xr(1,1) <= grid_1d.xmax)
    if (xr(1,1) >= grid_1d.xmin) && (xr(1,2) <= grid_1d.xmax)  
he(ifig) =  plot(xr(1,:),[pstates(ifig,1),pstates(ifig,1)],'k','LineWidth',LW);  
    elseif (xr(1,1) <= grid_1d.xmin) && (xr(1,2) <= grid_1d.xmax)      
he(ifig) =  plot([grid_1d.xmin,xr(1,2)],[pstates(ifig,1),pstates(ifig,1)],'k','LineWidth',LW);    
    elseif (xr(1,1) >= grid_1d.xmin) && (xr(1,2) >= grid_1d.xmax)  
he(ifig) =  plot([xr(1,1),grid_1d.xmax],[pstates(ifig,1),pstates(ifig,1)],'k','LineWidth',LW);    
    end
    end
hold on
  for i=2:nregio
    if (xr(i,2) >= grid_1d.xmin) && (xr(i,1) <= grid_1d.xmax)
    if (xr(i,1) >= grid_1d.xmin) && (xr(i,2) <= grid_1d.xmax)  
he(ifig) =  plot(xr(i,:),[pstates(ifig,i),pstates(ifig,i)],'k','LineWidth',LW);  
    elseif (xr(i,1) <= grid_1d.xmin) && (xr(i,2) <= grid_1d.xmax)      
he(ifig) =  plot([grid_1d.xmin,xr(i,2)],[pstates(ifig,i),pstates(ifig,i)],'k','LineWidth',LW);    
    elseif (xr(i,1) >= grid_1d.xmin) && (xr(i,2) >= grid_1d.xmax)  
he(ifig) =  plot([xr(i,1),grid_1d.xmax],[pstates(ifig,i),pstates(ifig,i)],'k','LineWidth',LW);    
    end
    end
  end  

  ifl = find(xfl >= grid_1d.xmin & xfl <= grid_1d.xmax);
  ial = find(xal >= grid_1d.xmin & xal <= grid_1d.xmax);
  isl = find(xsl >= grid_1d.xmin & xsl <= grid_1d.xmax);
  icd = find(xcd >= grid_1d.xmin & xcd <= grid_1d.xmax);
  isr = find(xsr >= grid_1d.xmin & xsr <= grid_1d.xmax);
  iar = find(xar >= grid_1d.xmin & xar <= grid_1d.xmax);
  ifr = find(xfr >= grid_1d.xmin & xfr <= grid_1d.xmax);

  plot(xfl(ifl),qfl(ifig,ifl),'k','LineWidth',LW);
  plot(xal(ial),qal(ifig,ial),'k','LineWidth',LW);
  plot(xsl(isl),qsl(ifig,isl),'k','LineWidth',LW);
  plot(xcd(icd),qcd(ifig,icd),'k','LineWidth',LW);
  plot(xsr(isr),qsr(ifig,isr),'k','LineWidth',LW);
  plot(xar(iar),qar(ifig,iar),'k','LineWidth',LW);
  plot(xfr(ifr),qfr(ifig,ifr),'k','LineWidth',LW);

elseif ifig == 8
  figure(ifig)
%*   clf
    if (xr(1,2) >= grid_1d.xmin) && (xr(1,1) <= grid_1d.xmax)
    if (xr(1,1) >= grid_1d.xmin) && (xr(1,2) <= grid_1d.xmax)  
he(ifig) =  plot(xr(1,:),[states(5,1),states(5,1)],'k','LineWidth',LW);  
    elseif (xr(1,1) <= grid_1d.xmin) && (xr(1,2) <= grid_1d.xmax)      
he(ifig) =  plot([grid_1d.xmin,xr(1,2)],[states(5,1),states(5,1)],'k','LineWidth',LW);    
    elseif (xr(1,1) >= grid_1d.xmin) && (xr(1,2) >= grid_1d.xmax)  
he(ifig) =  plot([xr(1,1),grid_1d.xmax],[states(5,1),states(5,1)],'k','LineWidth',LW);    
    end
    end

  hold on

  for i=2:nregio
    if (xr(i,2) >= grid_1d.xmin) && (xr(i,1) <= grid_1d.xmax)
    if (xr(i,1) >= grid_1d.xmin) && (xr(i,2) <= grid_1d.xmax)  
he(ifig) =  plot(xr(i,:),[states(5,i),states(5,i)],'k','LineWidth',LW);  
    elseif (xr(i,1) <= grid_1d.xmin) && (xr(i,2) <= grid_1d.xmax)      
he(ifig) =  plot([grid_1d.xmin,xr(i,2)],[states(5,i),states(5,i)],'k','LineWidth',LW);    
    elseif (xr(i,1) >= grid_1d.xmin) && (xr(i,2) >= grid_1d.xmax)  
he(ifig) =  plot([xr(i,1),grid_1d.xmax],[states(5,i),states(5,i)],'k','LineWidth',LW);    
    end
    end
  end  

  plot(xfl(ifl),ufl(5,ifl),'k','LineWidth',LW);
  plot(xal(ial),ual(5,ial),'k','LineWidth',LW);
  plot(xsl(isl),usl(5,isl),'k','LineWidth',LW);
  plot(xcd(icd),ucd(5,icd),'k','LineWidth',LW);
  plot(xsr(isr),usr(5,isr),'k','LineWidth',LW);
  plot(xar(iar),uar(5,iar),'k','LineWidth',LW);
  plot(xfr(ifr),ufr(5,ifr),'k','LineWidth',LW);
  
elseif ifig == 9
  figure(ifig)
%*   clf
    if (xr(1,2) >= grid_1d.xmin) && (xr(1,1) <= grid_1d.xmax)
    if (xr(1,1) >= grid_1d.xmin) && (xr(1,2) <= grid_1d.xmax)  
he(ifig) =  plot(xr(1,:),[psi(1),psi(1)],'k','LineWidth',LW);  
    elseif (xr(1,1) <= grid_1d.xmin) && (xr(1,2) <= grid_1d.xmax)      
he(ifig) =  plot([grid_1d.xmin,xr(1,2)],[psi(1),psi(1)],'k','LineWidth',LW);    
    elseif (xr(1,1) >= grid_1d.xmin) && (xr(1,2) >= grid_1d.xmax)  
he(ifig) =  plot([xr(1,1),grid_1d.xmax],[psi(1),psi(1)],'k','LineWidth',LW);    
    end
    end

  hold on

  for i=2:nregio
    if (xr(i,2) >= grid_1d.xmin) && (xr(i,1) <= grid_1d.xmax)
    if (xr(i,1) >= grid_1d.xmin) && (xr(i,2) <= grid_1d.xmax)  
he(ifig) =  plot(xr(i,:),[psi(i),psi(i)],'k','LineWidth',LW);  
    elseif (xr(i,1) <= grid_1d.xmin) && (xr(i,2) <= grid_1d.xmax)      
he(ifig) =  plot([grid_1d.xmin,xr(i,2)],[psi(i),psi(i)],'k','LineWidth',LW);    
    elseif (xr(i,1) >= grid_1d.xmin) && (xr(i,2) >= grid_1d.xmax)  
he(ifig) =  plot([xr(i,1),grid_1d.xmax],[psi(i),psi(i)],'k','LineWidth',LW);    
    end
    end
  end  
  
  psifl = atan2(qfl(7,:),qfl(6,:));
  psial = atan2(qal(7,:),qal(6,:));
  psisl = atan2(qsl(7,:),qsl(6,:));
  psicd = atan2(qcd(7,:),qcd(6,:));
  psisr = atan2(qsr(7,:),qsr(6,:));
  psiar = atan2(qar(7,:),qar(6,:));
  psifr = atan2(qfr(7,:),qfr(6,:));
  plot(xfl(ifl),psifl(ifl),'k','LineWidth',LW);
  plot(xal(ial),psial(ial),'k','LineWidth',LW);
  plot(xsl(isl),psisl(isl),'k','LineWidth',LW);
  plot(xcd(icd),psicd(icd),'k','LineWidth',LW);
  plot(xsr(isr),psisr(isr),'k','LineWidth',LW);
  plot(xar(iar),psiar(iar),'k','LineWidth',LW);
  plot(xfr(ifr),psifr(ifr),'k','LineWidth',LW);
end  


if ~strcmp(athena_output,'latex')
  ylim(limits_y(ifig,:));  
  ylabel(labels_y{ifig})
  title(['Problem: ',prob])
else
%*   if (ifig == 1)
%*     text(0.02,0.9,labels_y{ifig})
%*   else
%*     text(0.02,0.6,labels_y{ifig})  
%*   end

  if ~exist('set_ylim','var')
    set_ylim = 'off';
  end

  if strcmp(set_ylim,'on')
    ylim(limits_y(ifig,:));
    if strcmp(plot1,'plot_dissertation')
      ylabel(labels_y{ifig})
    end
  end

  xlim([xstart xend]);  
end