%******************************************************************************* 
%* Program: output_latex_table.m
%* Description: Out solution to latex table format
%* Author: Andrew Kercher 
%* References: 
%******************************************************************************* 

%-------------------------------------------------------------------------------
% Create table with solution
%-------------------------------------------------------------------------------
fname = sprintf('%s/dissertation/tab/%s_table.tex',pwd,prob);
if(strcmp(prob,'im24'))
  fname = sprintf('%s/dissertation/tab/%s_table.tex',pwd,test);  
elseif(strcmp(prob,'AK6'))
  fname = sprintf('%s/dissertation/tab/%s_table.tex',pwd,test);  
end  

fid = fopen(fname,'w');


fprintf(fid,'\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}} ccccccc}\n');
fprintf(fid,'\\\\ \n');
fprintf(fid,'\\hline \n'); 
fprintf(fid,'\\hline \n');
%* fprintf(fid,'$\\rho$ & $v_n$ & $v_y$ & $v_z$ & $p_g$ & $B_y$ & $B_z$ \\\\ \n');
fprintf(fid,'$\\rho$ & $v_n$ & $v_y$ & $v_z$ & $p_g$ & $B_t$ & $\\psi$ \\\\ \n');
fprintf(fid,'\\hline \n');

a = zeros(nunks,1);
b = zeros(nunks,1);
for j=1:nregio
  a(:) = 0.0;
  b(:) = 0;
  for i = 1:nunks
    a(i) = pstates(i,j);
    if i == 6
      a(i) = bperp(j);
    elseif i == 7
      a(i) = psi(j);  
    end

    if (abs(a(i)) > 1.0e-9)
%*     if abs(a(i)) < 1.0
    if (1.0 - abs(a(i))) > 1.0e-6  
      while ((abs(a(i)) < 1.0) && (b(i) > -20))
        b(i) -= 1;
        a(i) *= 10;
      end
    elseif abs(a(i)) >= 10.0
        while ((abs(a(i)) < 1.0) && (b(i) < 100))
        b(i) += 1;
        a(i) /= 10;
	end    
    end 
    else
	a(i) = 0.0;	
    end
  end
  
  ac = cellstr(num2str(a));
  
  for i = 1:nunks
    if a(i) >= 0  
      ac(i) = sprintf('\\phantom{-}% -0.4f\\mbox{\\sc{e}\\smaller{+}}%d',a(i),abs(b(i)));  
      if b(i) < 0
        ac(i) = sprintf('\\phantom{-}% -0.4f\\mbox{\\sc{e}\\smaller{-}}%d',a(i),abs(b(i)));        
      end
    else
      ac(i) = sprintf('% -0.4f\\mbox{\\sc{e}\\smaller{+}}%d',a(i),abs(b(i)));  
      if b(i) < 0
        ac(i) = sprintf('% -0.4f\\mbox{\\sc{e}\\smaller{-}}%d',a(i),abs(b(i)));          
      end
    end
  end

  fprintf(fid,'$%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ \\\\ \n',...
            ac{1},ac{2},ac{3},ac{4},ac{5},ac{6},ac{7});
  
end

fprintf(fid,'\\hline \n');
fprintf(fid,'\\end{tabular*} \n');

fclose(fid);

%-------------------------------------------------------------------------------
% Create table with wave speeds and error
%-------------------------------------------------------------------------------

fname = sprintf('%s/dissertation/tab/%s_err_table.tex',pwd,prob);
if(strcmp(prob,'im24'))
  fname = sprintf('%s/dissertation/tab/%s_err_table.tex',pwd,test);  
elseif(strcmp(prob,'AK6'))
  fname = sprintf('%s/dissertation/tab/%s_err_table.tex',pwd,test);  
end  
fid = fopen(fname,'w');

fprintf(fid,'\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}} cccccc}\n');
fprintf(fid,'\\\\ \n');
fprintf(fid,'\\hline \n'); 
fprintf(fid,'\\hline \n');
fprintf(fid,'Iteration & Residual & $W_{fl}$ & $W_{sl}$ & $W_{sr}$ & $W_{fr}$ \\\\ \n');
fprintf(fid,'\\hline \n');


iend = find(latex_ws(:,2) < 1.0e-6);
jend = iend(1);

if jend > 5
  jstart = jend - 4;
else
  jstart = 1;
  jend = 5;
end

a = zeros(6,1);
b = zeros(6,1);
for j=jstart:jend
  a(:) = 0.0;
  b(:) = 0;
  for i = 1:6
    a(i) = latex_ws(j,i);
    if (abs(a(i)) > 1.0e-12)
    if abs(a(i)) < 1.0
      while ((abs(a(i)) < 1.0) && (b(i) > -20))
        b(i) -= 1;
        a(i) *= 10;
      end
    elseif abs(a(i)) > 10.0
        while ((abs(a(i)) < 1.0) && (b(i) < 100))
        b(i) += 1;
        a(i) /= 10;
      end
    end 
  end
  end
  
  ac = cellstr(num2str(a));
  
  for i = 1:6
    if a(i) >= 0  
      ac(i) = sprintf('\\phantom{-}% -0.4f\\mbox{\\sc{e}\\smaller{+}}%d',a(i),abs(b(i)));  
      if b(i) < 0
        ac(i) = sprintf('\\phantom{-}% -0.4f\\mbox{\\sc{e}\\smaller{-}}%d',a(i),abs(b(i)));        
      end
    else
      ac(i) = sprintf('% -0.4f\\mbox{\\sc{e}\\smaller{+}}%d',a(i),abs(b(i)));  
      if b(i) < 0
        ac(i) = sprintf('% -0.4f\\mbox{\\sc{e}\\smaller{-}}%d',a(i),abs(b(i)));          
      end
    end
  end

  it_num = a(1) - jstart + 1;
  fprintf(fid,'$%d$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ \\\\ \n',...
            it_num,ac{2},ac{3},ac{4},ac{5},ac{6});
  
end

fprintf(fid,'\\hline \n');
fprintf(fid,'\\end{tabular*} \n');

fclose(fid);
