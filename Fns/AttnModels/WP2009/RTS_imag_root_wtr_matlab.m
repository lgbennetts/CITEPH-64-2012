function [x,exitflag]=RTS_imag_root_wtr_matlab(lam,H,i0,i1)

%% convergenc parameters:
EPS=1e-9;
MAXIT=100;

fl=lam*i0*sin(i0)+H*cos(i0);
fh=lam*i1*sin(i1)+H*cos(i1);

%% check there is a root in interval
if (fl > 0 & fh > 0) | (fl < 0 & fh < 0 )
  x=0;
  exitflag=0;
  return;
end

%% if either endpoint is "close enough", just use that for root
if (fl == 0)
  x = i0;
  exitflag=1;
  return;
end
if (fh == 0)
  x = i1;
  exitflag=2;
  return;
end

%% define search direction to be in dirn of increasing p
if (fl < 0.0)
  xl=i0;
  xh=i1;
else
  xl=i1;
  xh=i0;
end

%% initial guess for NR:
x=0.5*(i0+i1);
dxold=i1-i0;
dx=dxold;
p=lam*x*sin(x)+H*cos(x);
dp=(lam-H)*sin(x)+lam*x*cos(x);

%% find a root using N-R:*/
for its=1:MAXIT
  crit1=( ((x-xh)*dp-p)*((x-xl)*dp-p) >= 0 );
  crit2=( abs(2*p) > abs(dxold*dp) );
  if crit1 | crit2
    %% bisect if NR jumps out of range,
    %% or if not decreasing fast enough.*/
    dxold=dx;
    dx=0.5*(xh-xl);
    x=xl+dx;
    if (x == xl)
      %% change in root is negligible.
      exitflag=3;
      return;
    end
  else
    %% use NR
    dxold=dx;
    dx=p/dp;
    x0=x;
    x=x0-dx;
    if x == x0
      exitflag=3;
      return;
      %% change in root is negligible.
    end
  end
  if (abs(dx) < EPS)
    exitflag=3;
    return;
  end
  p=lam*x*sin(x)+H*cos(x);
  dp=(lam-H)*sin(x)+lam*x*cos(x);
  if (p < 0.0 )%% maintain the bracket on the root
    xl=x;
  else
    xh=x;
  end
end
x=0;
exitflag=4;
return;
%% if program reaches here it
%% has done MAXIT runs but hasn't converged yet