%******************************************************************************* 
%* Program: cons2prim_mhd_1d.m
%* Description: Converts conservative variables to primitive variables.
%* Author: Andrew Kercher
%-------------------------------------------------------------------------------
%* Input:
%*   nunks: number of unknowns
%*   ndof : number of degrees of freedom
%*   gamma: gas constant
%*   qcons: conservative variables (d,d*vx,d*vy,d*vz,en,by,bz)
%*   bx   : x-component of magnetic field  
%* Output:
%*   qprim: primitive variables (d,vx,vy,vz,pg,by,bz)
%******************************************************************************* 


function qprim = cons2prim_mhd_1d(nunks,ndof,gamma,qcons,bx)
  
%*   qprim = zeros(nunks,ndof);  
  
  for idof=1:ndof    
    vx  = qcons(2,idof)/qcons(1,idof);
    vy  = qcons(3,idof)/qcons(1,idof);    
    vz  = qcons(4,idof)/qcons(1,idof);          
    v2  = vx*vx + vy*vy + vz*vz;
    by  = qcons(6,idof);
    bz  = qcons(7,idof);      
    b2  = bx*bx + by*by + bz*bz;
  
    pg = max(0,(gamma - 1)*(qcons(5,idof) - 0.5*(qcons(1,idof)*v2 + b2)));
  
    qprim(1,idof) = qcons(1,idof);
    qprim(2,idof) = vx;
    qprim(3,idof) = vy;
    qprim(4,idof) = vz;      
    qprim(5,idof) = pg;
    qprim(6,idof) = qcons(6,idof);
    qprim(7,idof) = qcons(7,idof);    
      
  end