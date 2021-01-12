% Grafica la función de Rosenbrock:
%
%  f(x) = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2
%
% con único mínimo en x* = [1, 1]'.
% La graficación es en la caja [0.8, 1.2]x[0.8, 1.2].
%
% Análisis Aplicado
% ITAM
%-------------------------------------------------------
close all
fname = 'rosenbrock';
x1 = 0.1; x2 = 1.9;
y1 = 0.1; y2 = 1.9;
n = 50;
xx = linspace(x1,x2,n)';
yy = linspace(y1,y2,n)';
Z = zeros(n);

for i = 1:n
    a = xx(i);
    for j = 1:n
        Z(i,j) = feval(fname,[a,yy(j)]');
    end
end

[XX,YY]= meshgrid(xx,yy);


mesh(XX,YY,Z)
title('Rosenbrock' ,'Fontsize',20)
xlabel('EJE X1','Fontsize',16)
ylabel('EJE X2','Fontsize',16)
hold on
plot3(1,1,0,'dr','linewidth',6)

