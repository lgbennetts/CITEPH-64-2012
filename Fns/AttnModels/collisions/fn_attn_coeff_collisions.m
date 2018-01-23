%% fn_attn_coeff_collisions.m
%% Author: Timothy Williams
%% Date: 20140910, 15:40:45 CEST
%%
%% out = attn_collisions(...
%%          wave_pram,ice_pram,RAOsurge,conc,sep)
%%
%% NB attenuation dependent on wave amplitude;
%% wave_pram = [om,A,k,cg]: radial freq, amplitude, wave number, group vel
%% ice_pram  = [D,h,rhoi]:  ice floe diameter, thickness and density
%% RAOsurge: surge response for 1m amplitude wave
%% conc: ice concentration
%% sep: ice floe separation (distance between centres)
%% NB sep>D.
%% TODO: calc sep from conc & D if it's not specified
%%
%% out = struct:
%%   out.atten_coeff - atten coeff
%%   out.status      - any collisions (0/1)
%%   out.omt1        - time of collision 1
%%   out.omt2        - time of collision 2
%%   out.V1          - velocity of collision 1
%%   out.V2          - velocity of collision 2

function out = fn_attn_coeff_collisions(...
   coll_inputs,wave_pram,ice_pram,RAOsurge,conc,sep,f_coll)

do_test  = 0;
g        = 9.81;
rhow     = 1025;

if nargin==0
   do_test  = 1;
   %%
   om          = 2*pi/14;
   k           = om^2/g;
   cp          = om/k;
   cg          = g/2/om;
   wave_pram   = [om,k,cg];
   %%
   D           = 99;
   sep         = 100;
   conc        = (pi*D^2/4)/(sep^2);
   h           = 1;
   rhoi        = 922.5;
   ice_pram    = [D,h,rhoi];
   %%
   A           = 1;
   rest_coeff  = .9;
   coll_inputs = [A,rest_coeff];
   %%
   RAOsurge = 1;
end


A          = coll_inputs(:,1);%%wave amplitudes
rest_coeff = coll_inputs(:,2);%%restitution coefficients
alp        = 0*A;%%Attenuation coefficient if no collision;
Ncoll      = length(A);
%%
omt1     = NaN+0*A;
omt2     = NaN+0*A;
V1       = NaN+0*A;
V2       = NaN+0*A;
status   = 0*A;
%%
om = wave_pram(1);%%wave radial frequency
k  = wave_pram(2);%%wave number
cg = wave_pram(3);%%group velocity

f  = om/2/pi;%%wave frequency
if ~exist('f_coll','var')
   f_coll   = 2*f;%%collision freq = 2x wave frequency
end

D     = ice_pram(1);%%floe diameter
h     = ice_pram(2);%%ice thickness
rhoi  = ice_pram(3);%%ice density

A_floe    = pi*D^2/4;
rho_floes = conc/A_floe;%%floes per unit area - only needed for testing, cancels out

if ~exist('conc','var')
   conc  = 1;
end

err.identifier = 'AttnModels:fn_atten_coeff_collisions';
if ~exist('sep','var')
   %% WORK OUT sep from conc & D
   %% TODO: implement this
   err.message = 'argument "sep" not given';
   error(err);
end
if sep<=D
   disp(['sep = ',num2str(sep),'m; D = ',num2str(D),'m']);
   err.message = ['separation of floe centres ("sep") ',...
                  'should be larger than the diameter ("D")'];
   error(err);
end

X    = RAOsurge*A;%%absolute surge
ks   = k*sep;
Csq  = 2*(1-cos(ks));
C    = sqrt(Csq);
bet0 = (rhoi*h*om^2)/(2*rhow*g)*RAOsurge^2;

verbose = 0;
if C==0
   %% independent of wave amplitude
   %% (only a function of wave freq and separation)
   if verbose
      disp('No collision - floes moving exactly in phase');
   end
   out.atten_coeff = alp;
   out.status = status;
   out.omt1   = omt1;
   out.omt2   = omt2;
   out.V1     = V1;
   out.V2     = V2;
   return
end

for jc=1:Ncoll

   %% check if floes moving exactly in phase
   %% (no collisions)

   %% check if surge amplitudes too small
   %% (no collisions)
   cosZ  = (sep-D)/(X(jc)*C);%%Z=om*t+phi
   if abs(cosZ)>1
      if verbose
         disp([jc,Ncoll]);
         disp(['No collision - surge amplitudes ',...
               'too small to produce a collision']);
      end
      continue;
   end

   %% collisions if got to here
   status(jc) = 1;

   %% angle phi - 2 possibilities
   %% NB phi2 = phi+pi gives the same dE
   tan_phi = 2*sin(ks)/Csq;
   phi     = atan(tan_phi);
   zet1    = acos(cosZ);
   
   bet     = bet0*(1-rest_coeff(jc)^2);%%nondimensional
   F_av    = fn_Fcoll(2*ks,2*zet1-2*phi);
   alp(jc) = (conc*f_coll*bet.*F_av)/cg;%%attenuation per meter

   %%extra outputs:
   %% times and speeds of collisions
   %% - for use in drag model
   omt1(jc) = zet1-phi;
   omt2(jc) = -zet1-phi;
   V1(jc)   = -om*X(jc)*sin(omt1(jc));
   V2(jc)   = -om*X(jc)*sin(omt2(jc));
end

out.atten_coeff = alp;
out.status = status;
out.omt1   = omt1;
out.omt2   = omt2;
out.V1     = V1;
out.V2     = V2;
if 0
   A,omt1,V1,omt2,V2
   pause
end

if do_test
   np2   = 100;
   omt   = linspace(-pi,pi,np2).';
   x1    = D/2+X(1)*cos(omt);
   x2    = sep-D/2+cos(ks-omt);
   V1    = -om*X(1)*sin(omt);
   V2    = -om*X(1)*sin(ks-omt);
   %%
   omt1  = zet1(1)-phi;
   x1_1  = D/2+X(1)*cos(omt1);
   x2_1  = sep-D/2+cos(ks-omt1);
   V1_1  = -om*X(1)*sin(omt1);
   V2_1  = -om*X(1)*sin(ks-omt1);
   [omt1/pi,V1_1^2+V2_1^2]
   %%
   omt2  = -zet1(1)-phi;
   x1_2  = D/2+X(1)*cos(omt2);
   x2_2  = sep-D/2+cos(ks-omt2);
   V1_2  = -om*X(1)*sin(omt2);
   V2_2  = -om*X(1)*sin(ks-omt2);
   [omt2/pi,V1_2^2+V2_2^2]
   %%
   subplot(2,1,1);
   plot(omt/pi,[x1,x2]);
   GEN_proc_fig('\omega{t}/\pi','x, m')
   hold on;
   plot(omt1/pi,[x1_1,x2_1],'^');
   plot(omt2/pi,[x1_2,x2_2],'o');
   hold off;
   %%
   subplot(2,1,2);
   plot(omt/pi,[V1,V2]);
   GEN_proc_fig('\omega{t}/\pi','V, ms^{-1}')
   hold on;
   plot(omt1/pi,[V1_1,V2_1],'^');
   plot(omt2/pi,[V1_2,V2_2],'o');
   hold off;

   %%test at omt1:
   rc    = rest_coeff(1);
   Vsq   = V1_1^2+V2_1^2;
   Usq   = rc^2*Vsq;
   m     = rhoi*h*A_floe;
   dE1   = .5*m*(Vsq-Usq);

   V  = om*RAOsurge*A(1);%X,RAOsurge*A
   E  = rhow*g*A(1)^2/2;
   %[Vsq,V^2*(sin(omt1)^2+sin(k*sep-omt1)^2),...
   % V^2/2*(2-cos(2*omt1)-cos(2*(k*sep-omt1))),V^2/2*F_av]

   %dE1
   %.5*A_floe*[rhoi*h*(1-rc^2)*Vsq,rhoi*h*(1-rc^2)*V^2/2*F_av,...
   %                            (rhoi*h*om^2)/(rhow*g)*(1-rc^2)*RAOsurge^2*F_av*E ]

   %%test at omt2:
   Vsq = V1_2^2+V2_2^2;
   Usq = rc^2*Vsq;
   dE2 = .5*m*(Vsq-Usq);
   %%
   if 0
      dE          = [dE1,dE2]
      dE_theory   = bet*A_floe*E*F_av

      %%
      Phi   = f_coll*rho_floes*[dE1,dE2];
      %f_coll*conc*bet*E*F_av,alp*cg*E
      (Phi/E/cg),alp
   end

   if 0
      F1 = 2-cos(2*omt1)-cos(2*omt1-2*ks);
      F2 = 2-cos(2*omt2)-cos(2*omt2-2*ks);
      [F1,F2,F_av]
      %%
      sin(2*zet1(1))*(sin(2*phi)+sin(2*(phi+ks)))
      sin(2*zet1(1)),sin(2*phi)+sin(2*(phi+ks))
      4*tan(phi)*cos(ks)+2*sin(ks)*(1-tan(phi)^2)
      %%
      s  = sin(ks);
      c  = cos(ks);
      4*s*c*(1-c)+2*s*(1-c)^2-2*s*s^2
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F=fn_Fcoll(xi,zet)

%F  = (3-cos(xi))*(1-cos(zet))+...
%     (1-cos(xi))*(1+cos(zet))+...
%     -2*sin(xi)*sin(zet);
F  = 2-cos(zet)-cos(zet-xi);
