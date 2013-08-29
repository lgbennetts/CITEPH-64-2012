% function tvec = fn_TestTimes(f)
%
% INPUTS:
%
% f = wave frequency [Hz]
% x = location of point in question
%
% OUTPUT:
%
% Tindic = structure containing indicitive time
%          field are 'time' & 'description'
%
% written by L Bennetts 01.08.2013 / Oceanide, La Seyne

function Tindic = fn_TestTimes(f,x)

% Waves

Param = ParamDef3d_Oceanide;

fortyp='freq'; 
bed=3;

Forcing = Force_def(Param.g, bed, fortyp, f);

cg = fn_grpvel(2*pi/Forcing.lam0,bed,Forcing.f);
cp = f*Forcing.lam0;

% Geometry [m] (origin is middle of tank)

WMK = -16;

MIZ = [-2.5,2.5];
 
BCH = 16;

% Significant times [s]

t0 = 50;                   d0='wave maker switched on';
t1 = t0 + (MIZ(1)-WMK)/cg; d1='waves reach MIZ';
t2 = t0 + (MIZ(2)-WMK)/cg; d2='wave exit MIZ';
t3 = t0 + (BCH-WMK)/cg;    d3='wave reach beach';
t4 = t0 + 2*(BCH-WMK)/cg;  d4='reflected wave from beach reach wave maker';

% times relative to specified point [s]

if x<MIZ(1) % point on LHS

 tx0 = t0+(x-WMK)/cg;                   dx0='waves reach x';
 tx1 = t0+((MIZ(1)-x)+(MIZ(1)-WMK))/cg; dx1='waves ref by MIZ reach x';
 tx2 = t0+((BCH-x)+(BCH-WMK))/cg;       dx2='waves ref beach MIZ reach x';

elseif x>MIZ(2) % point on RHS

 tx0 = t0+(x-WMK)/cg;                   dx0='waves reach x';
 tx1 = t0+((BCH-x)+(BCH-WMK))/cg;       dx1='waves reflected by beach reach x';
 tx2 = tx1+2*(x-MIZ(2));                dx2='waves re-ref by MIZ reach x';
 
end

t = {t0,t1,t2,t3,t4,tx0,tx1,tx2};
des = {d0,d1,d2,d3,d4,dx0,dx1,dx2};


Tindic = struct('description',des,...
 'time',t);

return

% Group velocity

function cg = fn_grpvel(k,h,f)

% see Mei book p.18

cg = 0.5*(2*pi*f/k)*(1+(2*k*h/sinh(2*k*h)));

return