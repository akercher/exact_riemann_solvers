%******************************************************************************* 
%* Program: hlld.m
%* Description: Computes states and fluxes using HLLD approximation.
%*              This is an octave version of the "hlld.c" function
%*              included with Athena, an open source approximate MHD
%*              solver [2].  The source code is available for download
%*              via the projects homepage:
%*              https://trac.princeton.edu/Athena/
%* Author: Andrew Kercher 
%* References: 
%*         [1] T. Miyoshi & K. Kusano, "A multi-state HLL approximate
%*             Riemann solver for ideal MHD", JCP, 208, 315 (2005).
%*         [2] J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon
%*             "Athena: A new code for astrophysical MHD", ApJS, (2008).
%-------------------------------------------------------------------------------
%******************************************************************************* 

function [qhlld,fhlld,vspd] = hlld(nunks,gamma,bx,ql,qr);

spd = zeros(5,1);
qhlld = zeros(7,6);

% initialize intermediate states in primitive variables
qls = zeros(7,1);
qlss = zeros(7,1);
qrss = zeros(7,1);
qrs = zeros(7,1);

% initialize states in conservative variables
ul = zeros(7,1);
uls = zeros(7,1);
ulss = zeros(7,1);
urss = zeros(7,1);
urs = zeros(7,1);
ur = zeros(7,1);

ul(:) = prim2cons_mhd_1d(nunks,1,gamma,ql,bx);
ur(:) = prim2cons_mhd_1d(nunks,1,gamma,qr,bx);

% initialize fluxes
fl = zeros(7,1);
fr = zeros(7,1);
fhlld = zeros(7,1);

dl = ql(1);
sqrtdl = sqrt(dl);
Vl = 1/dl;
vnl = ql(2);
vyl = ql(3);
vzl = ql(4);
pgl = ql(5);
byl = ql(6);
bzl = ql(7);
btl = sqrt(ql(6)*ql(6) + ql(7)*ql(7));
kel = 0.5*(vnl*vnl + vyl*vyl + vzl*vzl);
pbl = 0.5*(bx*bx + btl*btl);
enl = pgl/(gamma - 1) + dl*kel + pbl;
ptl = pgl + pbl;
psil = atan2(bzl,byl);

dr = qr(1);
sqrtdr = sqrt(dr);
Vr = 1/dr;
vnr = qr(2);
vyr = qr(3);
vzr = qr(4);
pgr = qr(5,end);
byr = qr(6);
bzr = qr(7);
btr = sqrt(qr(6)*qr(6) + qr(7)*qr(7));
ker = 0.5*(vnr*vnr + vyr*vyr + vzr*vzr);
pbr = 0.5*(bx*bx + btr*btr);
enr = pgr/(gamma - 1) + dr*ker + pbr;
ptr = pgr + pbr;
psir = atan2(bzr,byr);

% left fluxes
fl(1) = ul(2);
fl(2) = ul(2)*ql(2) + ptl - bx*bx;
fl(3) = ul(2)*ql(3) - bx*ql(6);
fl(4) = ul(2)*ql(4) - bx*ql(7);
fl(5) = ql(2)*(ul(5) + ptl - bx*bx) - bx*(ql(3)*ql(6) + ql(4)*ql(7));
fl(6) = ql(2)*ql(6) - bx*ql(3);
fl(7) = ql(2)*ql(7) - bx*ql(4);

% right fluxes
fr(1) = ur(2);
fr(2) = ur(2)*qr(2) + ptr - bx*bx;
fr(3) = ur(2)*qr(3) - bx*qr(6);
fr(4) = ur(2)*qr(4) - bx*qr(7);
fr(5) = qr(2)*(ur(5) + ptr - bx*bx) - bx*(qr(3)*qr(6) + qr(4)*qr(7));
fr(6) = qr(2)*qr(6) - bx*qr(3);
fr(7) = qr(2)*qr(7) - bx*qr(4);

% Sound, Alfven, fast and slow characteristics speeds in Lagrangian mass coordiantes
[C0l2,Cal2,Ctl2,Csl2,Cfl2] = Lagrangian_wave_speeds_mhd(gamma,dl,pgl,btl,bx);
c0l = sqrt(C0l2)/dl;
cal = sqrt(Cal2)/dl;
ctl = sqrt(Ctl2)/dl;
csl = sqrt(Csl2)/dl;
cfl = sqrt(Cfl2)/dl;

% Sound, Alfven, fast and slow characteristics speeds in Lagrangian mass coordiantes
[C0r2,Car2,Ctr2,Csr2,Cfr2] = Lagrangian_wave_speeds_mhd(gamma,dr,pgr,btr,bx);
c0r = sqrt(C0r2)/dr;
car = sqrt(Car2)/dr;
ctr = sqrt(Ctr2)/dr;
csr = sqrt(Csr2)/dr;
cfr = sqrt(Cfr2)/dr;

cfmax = max(cfl,cfr);

if vnl <= vnr
  spd(1) = vnl - cfmax;
  spd(5) = vnr + cfmax;  
else
  spd(1) = vnr - cfmax;
  spd(5) = vnl + cfmax;    
end  

if spd(1) >= 0 
  fhlld = fl;
  return;
end

if spd(5) <= 0 
  fhlld = fr;
  return;
end

% middle speeds
sdl = spd(1) - vnl;
sdr = spd(5) - vnr;

spd(3) = (sdr*dr*vnr - sdl*dl*vnl - ptr + ptl)/(sdr*dr - sdl*dl);

sdml = spd(1) - spd(3);
sdmr = spd(5) - spd(3);

qls(1) = dl*sdl/sdml;
qrs(1) = dr*sdr/sdmr;

sqrtdls = sqrt(qls(1));
sqrtdrs = sqrt(qrs(1));

spd(2) = spd(3) - abs(bx)/sqrtdls;
spd(4) = spd(3) + abs(bx)/sqrtdrs;

% left intermediate state

pts = ptl + dl*sdl*(sdl-sdml);

qls(2) = spd(3);

tmp = bx*(sdl-sdml)/(dl*sdl*sdml - bx*bx);
qls(3) = vyl - byl*tmp;
qls(4) = vzl - bzl*tmp;

tmp = (dl*sdl*sdl - bx*bx)/(dl*sdl*sdml - bx*bx);
qls(6) = byl*tmp;
qls(7) = bzl*tmp;

vbl = vnl*bx + vyl*byl + vzl*bzl;
vbls = qls(2)*bx + qls(3)*qls(6) + qls(4)*qls(7);

enls = (sdl*ul(5) - ptl*vnl + pts*spd(3) ...
       + bx*(vbl - vbls))/sdml;

uls(1) = qls(1);
uls(2) = qls(1)*qls(2);
uls(3) = qls(1)*qls(3);
uls(4) = qls(1)*qls(4);
uls(5) = enls;
uls(6) = qls(6);
uls(7) = qls(7);

qls = cons2prim_mhd_1d(nunks,1,gamma,uls,bx);

% right intermediate state
qrs(2) = spd(3);

tmp = bx*(sdr-sdmr)/(dr*sdr*sdmr - bx*bx);
qrs(3) = vyr - byr*tmp;
qrs(4) = vzr - bzr*tmp;

tmp = (dr*sdr*sdr - bx*bx)/(dr*sdr*sdmr - bx*bx);
qrs(6) = byr*tmp;
qrs(7) = bzr*tmp;

vbr = vnr*bx + vyr*byr + vzr*bzr;
vbrs = qrs(2)*bx + qrs(3)*qrs(6) + qrs(4)*qrs(7);

enrs = (sdr*ur(5) - ptr*vnr + pts*spd(3) ...
       + bx*(vbr - vbrs))/sdmr;

urs(1) = qrs(1);
urs(2) = qrs(1)*qrs(2);
urs(3) = qrs(1)*qrs(3);
urs(4) = qrs(1)*qrs(4);
urs(5) = enrs;
urs(6) = qrs(6);
urs(7) = qrs(7);

qrs = cons2prim_mhd_1d(nunks,1,gamma,urs,bx);

% Middle states
sgnbx = bx/abs(bx);

invd = 1/(sqrtdl + sqrtdr);

qlss(1) = qls(1);
qlss(2) = qls(2);

qrss(1) = qrs(1);
qrss(2) = qrs(2);

qlss(3) = invd*(sqrtdl*qls(3) + sqrtdr*qrs(3) + sgnbx*(qrs(6) - qls(6)));
qrss(3) = qlss(3);

qlss(4) = invd*(sqrtdl*qls(4) + sqrtdr*qrs(4) + sgnbx*(qrs(7) - qls(7)));
qrss(4) = qlss(4);

qlss(6) = invd*(sqrtdl*qrs(6) + sqrtdr*qls(6) ... 
		+ sgnbx*sqrtdl*sqrtdr*(qrs(3) - qls(3)));
qrss(6) = qlss(6);

qlss(7) = invd*(sqrtdl*qrs(7) + sqrtdr*qls(7) ... 
		+ sgnbx*sqrtdl*sqrtdr*(qrs(4) - qls(4)));
qrss(7) = qlss(7);

vbss = spd(3)*bx + qlss(3)*qlss(6) + qlss(4)*qlss(7);
ulss(5) = uls(5) - sqrtdl*sgnbx*(vbls - vbss);
urss(5) = urs(5) + sqrtdr*sgnbx*(vbrs - vbss);

ulss(1) = qlss(1);
ulss(2) = qlss(1)*qlss(2);
ulss(3) = qlss(1)*qlss(3);
ulss(4) = qlss(1)*qlss(4);
ulss(6) = qlss(6);
ulss(7) = qlss(7);

qlss = cons2prim_mhd_1d(nunks,1,gamma,ulss,bx);

urss(1) = qrss(1);
urss(2) = qrss(1)*qrss(2);
urss(3) = qrss(1)*qrss(3);
urss(4) = qrss(1)*qrss(4);
urss(6) = qrss(6);
urss(7) = qrss(7);

qrss = cons2prim_mhd_1d(nunks,1,gamma,urss,bx);

qhlld(:,1) = ql;
qhlld(:,2) = qls;
qhlld(:,3) = qlss;
qhlld(:,end-2) = qrss;
qhlld(:,end-1) = qrs;
qhlld(:,end) = qr;


if spd(2) >= 0
  for i = 1:7
    fhlld(i) = fl(i) + spd(1)*(uls(i) - ul(i));
  end
elseif spd(3) >= 0
  for i = 1:7
    fhlld(i) = fl(i) - spd(1)*ul(i) - (spd(2) - spd(1))*uls(i) + spd(2)*ulss(i);  
  end
elseif spd(4) > 0
  for i = 1:7
    fhlld(i) = fr(i) - spd(5)*ur(i) - (spd(4) - spd(5))*urs(i) + spd(4)*urss(i);  
  end
else
  for i = 1:7
    fhlld(i) = fr(i) + spd(5)*(urs(i) - ur(i));
  end  
end

vspd = spd;
