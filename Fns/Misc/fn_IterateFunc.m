%Currently just iterated function for integers n

function out=fn_IterateFunc(func,n,x0,IPLOT)

if ~exist('func','var');    func= @sublinear; end
if ~exist('n','var');    n= 10; end
if ~exist('x0','var');   x0 = 1; end
if ~exist('IPLOT','var');   IPLOT = 1; end


close all;

Iterates = ones(n+1,1);
IteratesA = ones(n+1,1);
Iterates(1) = x0;
IteratesA(1) = x0;
for i  = 2:n+1
    Iterates(i) = func(Iterates(i-1));
    IteratesA(i) = linear(IteratesA(i-1));
end

if IPLOT == 1
   x = 0:x0/100.0:x0;
   yE = func(x);
   yA = linear(x);
   
   figure();
   p1 = plot(x,x,'-k', 'LineWidth',2,'DisplayName',' y = x');
   hold on;
   p2 = plot(x,yE,'-g', 'LineWidth',2,'DisplayName',' y = a*x - |f(x)|');
   p3 = plot(x,yA,'-b', 'LineWidth',2,'DisplayName',' y = a*x');
   
   PlotIteration(Iterates,'--r')
   p4 = plot(Iterates, func(Iterates),'.r','DisplayName','Real Iterate', 'MarkerSize',12);
   plot(func(Iterates), func(Iterates),'.r', 'MarkerSize',12);
   
   PlotIteration(IteratesA,':c')
   p5 = plot(IteratesA, linear(IteratesA),'.c','DisplayName','Approximate Iterate', 'MarkerSize',12);
   plot(linear(IteratesA), linear(IteratesA),'.c', 'MarkerSize',12);   
   
   xlabel('x')
   ylabel('y')
   title('Iteration Plot')
   legend([p1 p2 p3 p4 p5])
%    plot(x,y,'--k','DisplayName',' y = func(x)')
   
   
    
end


end

function PlotIteration(Iterates,LineSty)
    plot([Iterates(1),Iterates(1)], [0, Iterates(2)], LineSty)
for i = 1:length(Iterates)-2
    %Leg axis to function
    plot([Iterates(i),Iterates(i+1)], [Iterates(i+1), Iterates(i+1)], LineSty)
    plot([Iterates(i+1),Iterates(i+1)], [Iterates(i+1), Iterates(i+2)], LineSty)
end
end

function y = linear(x)
 y = 0.7*x;
end

function y = sublinear(x)
 y = 0.7*x - 0.1*x.^3;
end