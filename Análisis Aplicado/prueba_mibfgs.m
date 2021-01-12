%   Prueba_mibfgs.m
close all
f='Rosenbrock1000';
x0=(2.2)*ones(1000,1);
[x,k,Z]=mibfgs(f,x0);


v=[1:length(Z)]';
    semilogy(v,Z,'dr',v,Z,'b','Linewidth',3)
    xlabel('Iteraciones','Fontsize',16)
    ylabel('\|gx\|','Fontsize',16)
    title('Rosenbrock1000/BFGS/ Convergencia','Fontsize',16)