%%% generate some simple data to test WDM method
%%%
%%% written by L Bennetts Jan 2013 / Adelaide

% array of probes
A = [0 33:72:359]; dr=pi/180;  A = dr*A; % 28 Oct 87 DB  Corrected June 22 1994
R = [0 0.25 0.25 0.25 0.25 0.25] ;   % Radius of array (m). i.e. Polar coordinates of the staffs
X = R.*cos(A); Y = R.*sin(A);

ns=4; % sampling frequency
n=4096; % = 2^12

for run = 3
    
 if 0
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% simple harmonic wave
 
 th = pi/6; % angle wrt x-axis
 f = 1; % freq in Hz
 k0 = (2*pi*f)^2/9.81; % assume deep water
 
 phi = exp(1i*k0*(cos(th)*X+sin(th)*Y));
 
 t=0;
 for loop_t = 1:3*n
  data(loop_t,:) = real(phi*exp(-2i*pi*f*t));
  t=t+(1/ns);
 end
 
 description = ['simple harmonic wave: ang=' num2str(180*th/pi) '; f=' num2str(f)...
     ' Hz; k=' num2str(k0) ' 1/m'];
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 elseif 1
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% superposition of harmonic waves
 
 th = [-1,1]*pi/6; % angle wrt x-axis
 f = 1; % freq in Hz
 k0 = (2*pi*f)^2/9.81; % assume deep water
 
 phi = exp(1i*k0*(cos(th(1))*X+sin(th(1))*Y));
 for loop_a=2:length(th)
  phi = phi + exp(1i*k0*(cos(th(loop_a))*X+sin(th(loop_a))*Y));
 end
 
 t=0;
 for loop_t = 1:3*n
  data(loop_t,:) = real(phi*exp(-2i*pi*f*t));
  t=t+(1/ns);
 end
 
 angs = num2str(th(1));
 for loop_a=2:length(th); angs = [angs ', ' num2str(180*th(loop_a)/pi)]; end
         
 description = [int2str(length(th)) ' simple harmonic waves: angs=' ...
     angs '; f=' num2str(f) ' Hz; k=' num2str(k0) ' 1/m'];
 end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%% save data
    
 if run > 99
   eval(['save ../WDM/s13',int2str(run), ' data description'])
 elseif run > 9
  eval(['save ../WDM/s130',int2str(run), ' data description'])
 else
  eval(['save ../WDM/s1300',int2str(run), ' data description'])
 end
 
end