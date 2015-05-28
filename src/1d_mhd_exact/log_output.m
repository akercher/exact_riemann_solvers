%******************************************************************************* 
%* Program: log_output.m
%* Description: Writes to log file. 
%* Author: Andrew Kercher
%-------------------------------------------------------------------------------
%******************************************************************************* 

latex_ws = [latex_ws;kiter,err_k,wfl,wsl,wsr,wfr];

if abs(bx) > 0
  fprintf(flog,'\n');
  fprintf(flog,'iter = %d rfac_k = %04e err = %0.4e\n',kiter,rfac_k,err_k);  
  fprintf(flog,'wfl = % -0.4e wal = % -0.4e wsl = % -0.4e\n',wfl,wal,wsl);    
  fprintf(flog,'wfr = % -0.4e war = % -0.4e wsr = % -0.4e\n',wfr,war,wsr);      
  fprintf(flog,'dpsil = % -0.4e dpsil = % -0.4e\n',dpsil,dpsir);      
else
  fprintf(flog,'\n');
  fprintf(flog,'iter = %d rfac_k = %0.4e err = %0.4e\n',kiter,rfac_k,err_k);  
  fprintf(flog,'wmsl = % -0.4e wmsr = % -0.4e\n',wmsl,wmsr);    
end  

  fprintf(flog,'\n');
  fprintf(flog,'       d            vx            vy            vz            By            Bz            pg\n');
  fprintf(flog,' ------------------------------------------------------------------------------------------------\n');
  
  for i=1:nregio
  
    fprintf(flog,' % -0.4e   % -0.4e   % -0.4e   % -0.4e   % -0.4e   % -0.4e   % -0.4e\n',...
            pstates(1,i),pstates(2,i),pstates(3,i),pstates(4,i),...
            pstates(6,i),pstates(7,i),pstates(5,i));   
  
  end    
  fprintf(flog,' ------------------------------------------------------------------------------------------------\n\n');
