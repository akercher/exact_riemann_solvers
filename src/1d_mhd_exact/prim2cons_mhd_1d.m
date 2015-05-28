%******************************************************************************* 
%* Program: prim2cons_mhd_1d.m
%* Description: Converts primitive variables to conservative variables.
%* Author: Andrew Kercher
%-------------------------------------------------------------------------------
%* Input:
%*   nunks: number of unknowns
%*   ndof : number of degrees of freedom
%*   gamma: gas constant
%*   qprim: primitive variables (d,vx,vy,vz,pg,by,bz)
%*   bx   : x-component of magnetic field  
%* Output:
%*   qcons: conservative variables (d,d*vx,d*vy,d*vz,en,by,bz)
%******************************************************************************* 

function qcons = prim2cons_mhd_1d(nunks,ndof,gamma,qprim,bx)
  
%*   qcons = zeros(nunks,ndof);

  for idof=1:ndof
      vx = qprim(2,idof);
      vy = qprim(3,idof);      
      vz = qprim(4,idof);      
      v2 = vx*vx + vy*vy + vz*vz;
      by  = qprim(6,idof);
      bz  = qprim(7,idof);      

      b2  = bx*bx + by*by + bz*bz;
      
      qcons(1,idof) = qprim(1,idof);
      qcons(2,idof) = qprim(1,idof)*qprim(2,idof);    
      qcons(3,idof) = qprim(1,idof)*qprim(3,idof);        
      qcons(4,idof) = qprim(1,idof)*qprim(4,idof);              
      qcons(5,idof) = qprim(5,idof)/(gamma-1) + 0.5*(qprim(1,idof)*v2 + b2);
      qcons(6,idof) = qprim(6,idof);
      qcons(7,idof) = qprim(7,idof);

  end    
