% function tvec = fn_TestTimes(f)
%
% INPUTS:
%
% f   = wave frequency [Hz]
% x   = location of point in question
% TYP = 'attn' or 'single'
%
% OUTPUT:
%
% Tindic = structure containing indicitive time
%          field are 'time' & 'description'
%
% written by L Bennetts 01.08.2013 / Oceanide, La Seyne

function Tindic = fn_TestTimes(f,x,TYP)

if ~exist('f','var');   f=0.5; end
if ~exist('x','var');   x=-3; end
if ~exist('TYP','var'); TYP='single'; end

% Waves

Param = ParamDef_Oceanide;

fortyp='freq';

Forcing = Force_def(Param.g, Param.bed, fortyp, f);

cg = fn_grpvel(2*pi/Forcing.lam0,Param.bed,Forcing.f);
cp = f*Forcing.lam0;

%% ATTENUATION TESTS

if strcmp(TYP,'attn')
 
 % Geometry [m] (origin is middle of tank)
 
 WMK = -16;
 
 MIZ = [-2.5,2.5];
 
 BCH = 16;
 
 % Significant times [s]
 
 t0 = 50;                   d0='wave maker switched on'; c0='k:';
 t1 = t0 + (MIZ(1)-WMK)/cg; d1='waves reach MIZ'; c1='k:';
 t2 = t0 + (MIZ(2)-WMK)/cg; d2='waves leave MIZ'; c2='k:';
 t3 = t0 + (BCH-WMK)/cg;    d3='waves reach beach'; c3='k:';
 t4 = t0 + 2*(BCH-WMK)/cg;  d4='reflected waves from beach reach wave maker'; c4='k:';
 t5 = 120;                  d5='wave maker switched off'; c5='k:';
 
 % times relative to specified point [s]
 
 if x<MIZ(1) % point on LHS
  
  tx0 = t0+(x-WMK)/cg;                   dx0='waves reach x'; cx0='g:';
  tx1 = t0+((MIZ(1)-x)+(MIZ(1)-WMK))/cg; dx1='MIZ ref waves reach x'; cx1='b:';
  tx2 = t0+((BCH-x)+(BCH-WMK))/cg;       dx2='beach ref waves reach x'; cx2='m:';
  tx3 = t5+(x-WMK)/cg;                   dx3='final waves reach x'; cx3='r:';
  tx4 = t0+(2*(MIZ(1)-WMK)+(x-WMK))/cg;  dx4='wavemaker ref waves reach x'; cx4='c:';
  
 elseif x>MIZ(2) % point on RHS
  
  tx0 = t0+(x-WMK)/cg;                   dx0='waves reach x'; cx0='g:';
  tx1 = t0+((BCH-x)+(BCH-WMK))/cg;       dx1='beach ref waves reach x'; cx1='b:';
  tx2 = tx1+2*(x-MIZ(2));                dx2='MIZ re-ref waves reach x'; cx2='m:';
  tx3 = t5+(x-WMK)/cg;                   dx3='final waves reach x'; cx3='r:';
  
 else % within MIZ
  
  tx0 = t0+(x-WMK)/cg;                   dx0='waves reach x'; cx0='g:';
  tx1 = t0+((BCH-x)+(BCH-WMK))/cg;       dx1='beach ref waves reach x'; cx1='b:';
  tx2 = tx1+2*(x-MIZ(2));                dx2='MIZ re-ref waves reach x'; cx2='m:';
  tx3 = t5+(x-WMK)/cg;                   dx3='final waves reach x'; cx3='r:';
 
 end
 
 if exist('tx4','var')
  t = {t0,t1,t2,t3,t4,t5,tx0,tx1,tx2,tx3,tx4};
  des = {d0,d1,d2,d3,d4,d5,dx0,dx1,dx2,dx3,dx4};
  col = {c0,c1,c2,c3,c4,c5,cx0,cx1,cx2,cx3,cx4};
 else
  t = {t0,t1,t2,t3,t4,t5,tx0,tx1,tx2,tx3};
  des = {d0,d1,d2,d3,d4,d5,dx0,dx1,dx2,dx3};
  col = {c0,c1,c2,c3,c4,c5,cx0,cx1,cx2,cx3};
 end
 
 %% SINGLE FLOE TESTS
 
elseif strcmp(TYP,'single')
 
 % Geometry [m] (origin is middle of tank)
 
 WMK = -16;
 
 FLOE = [-0.495,.495];
 
 BCH = 16;
 
 % Significant times [s]
 
 t0 = 50;                    d0='wave maker switched on'; c0='k:';
 t1 = t0 + (FLOE(1)-WMK)/cg; d1='waves reach floe'; c1='k:';
 t2 = t0 + (FLOE(2)-WMK)/cg; d2='waves leave floe'; c2='k:';
 t3 = t0 + (BCH-WMK)/cg;     d3='waves reach beach'; c3='k:';
 t4 = t0 + 2*(BCH-WMK)/cg;   d4='reflected waves from beach reach wave maker'; c4='k:';
 t5 = 120;                   d5='wave maker switched off'; c5='k:';
 
 % times relative to specified point [s]
 
 if x<FLOE(1) % point on LHS
  
  tx0 = t0+(x-WMK)/cg;                     dx0='waves reach x';              cx0='g:';
  tx1 = t0+((FLOE(1)-x)+(FLOE(1)-WMK))/cg; dx1='waves ref by floe reach x';  cx1='b:';
  tx2 = t0+((BCH-x)+(BCH-WMK))/cg;         dx2='waves ref by beach reach x'; cx2='m:';
  tx3 = t5+(x-WMK)/cg;                     dx3='final waves reach x';        cx3='r:';
  
 elseif x>FLOE(2) % point on RHS
  
  tx0 = t0+(x-WMK)/cg;                    dx0='waves reach x';                    cx0='g:';
  tx2 = t0+((BCH-x)+(BCH-WMK))/cg;        dx2='waves reflected by beach reach x'; cx2='m:';
  tx1 = tx2+2*(x-FLOE(2));                dx1='waves re-ref by MIZ reach x';      cx1='b:';
  tx3 = t5+(x-WMK)/cg;                    dx3='final waves reach x';              cx3='r:';
  
 else % floe itself
  
  tx0 = t0+(FLOE(1)-WMK)/cg;              dx0='waves reach floe (x)';                         cx0='g:';
  tx1 = t0+3*(FLOE(1)-WMK)/cg;            dx1='waves reflected by wave maker reach floe (x)'; cx1='b:';
  tx2 = t0+((BCH-FLOE(2))+(BCH-WMK))/cg;  dx2='waves reflected by beach reach floe (x)';      cx2='m:';
  tx3 = t5+(FLOE(1)-WMK)/cg;              dx3='final waves reach floe (x)';                   cx3='r:';
  
 end

 t = {t0,t1,t2,t3,t4,t5,tx0,tx1,tx2,tx3};
 des = {d0,d1,d2,d3,d4,d5,dx0,dx1,dx2,dx3};
 col = {c0,c1,c2,c3,c4,c5,cx0,cx1,cx2,cx3};
 
 %% CALIBRATION TESTS
 
elseif strcmp(TYP,'calibration')
 
 % Geometry [m] (origin is middle of tank)
 
 WMK = -16;
 
 MIZ = [-2.5,2.5];
 
 BCH = 16;
 
 % Significant times [s]
 
 t0 = 50;                   d0='wave maker switched on'; c0='k:';
 t1 = t0 + (MIZ(1)-WMK)/cg; d1='waves reach MIZ'; c1='k:';
 t2 = t0 + (MIZ(2)-WMK)/cg; d2='waves leave MIZ'; c2='k:';
 t3 = t0 + (BCH-WMK)/cg;    d3='waves reach beach'; c3='k:';
 t4 = t0 + 2*(BCH-WMK)/cg;  d4='reflected waves from beach r each wave maker'; c4='k:';
 t5 = 120;                  d5='wave maker switched off'; c5='k:';

 
 if x<MIZ(1) % point on LHS
  
  tx0 = t0+(x-WMK)/cg;                   dx0='waves reach x'; cx0='g:';
  tx1 = t0+((MIZ(1)-x)+(MIZ(1)-WMK))/cg; dx1='MIZ ref waves reach x'; cx1='b:';
  tx2 = t0+((BCH-x)+(BCH-WMK))/cg;       dx2='beach ref waves reach x'; cx2='m:';
  tx3 = t5+(x-WMK)/cg;                   dx3='final waves reach x'; cx3='r:';
  
 elseif x>MIZ(2) % point on RHS
  
  tx0 = t0+(x-WMK)/cg;                   dx0='waves reach x'; cx0='g:';
  tx1 = t0+((BCH-x)+(BCH-WMK))/cg;       dx1='beach ref waves reach x'; cx1='b:';
  tx2 = tx1+2*(x-MIZ(2));                dx2='MIZ re-ref waves reach x'; cx2='m:';
  tx3 = t5+(x-WMK)/cg;                   dx3='final waves reach x'; cx3='r:';
  
 else % within MIZ
  
  tx0 = t0+(x-WMK)/cg;                   dx0='waves reach x'; cx0='g:';
  tx1 = t0+((BCH-x)+(BCH-WMK))/cg;       dx1='beach ref waves reach x'; cx1='b:';
  tx2 = tx1+2*(x-MIZ(2));                dx2='MIZ re-ref waves reach x'; cx2='m:';
  tx3 = t5+(x-WMK)/cg;                   dx3='final waves reach x'; cx3='r:';
 
 end
 
 t   = {t0,t1,t2,t3,t4,t5,tx0,tx1,tx2,tx3};
 des = {d0,d1,d2,d3,d4,d5,dx0,dx1,dx2,dx3};
 col = {c0,c1,c2,c3,c4,c5,cx0,cx1,cx2,cx3};
 
end

Tindic = struct('description',des,...
 'time',t,'colour',col);

return

% Group velocity

function cg = fn_grpvel(k,h,f)

% see Mei book p.18

cg = 0.5*(2*pi*f/k)*(1+(2*k*h/sinh(2*k*h)));

return