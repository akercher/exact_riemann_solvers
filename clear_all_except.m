%******************************************************************************* 
%* Program: clear_all_except
%* Description: clears all variables except ones specified.
%* Author: Andrew Kercher 
%-------------------------------------------------------------------------------
%******************************************************************************* 

save('./bin/save.tmp','var_list');
clear all
load('./bin/save.tmp');
com = sprintf('delete %s/bin/save.tmp',pwd);
eval(com);
