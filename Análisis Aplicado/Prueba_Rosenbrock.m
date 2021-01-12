% Prueba: Rosenbrock.m
close all
fname = 'rosenbrock';
x = [2,3]';

[x,iter] = MetodoBL(fname,x);
%[x,iter,W] = Metodo_Descenso_Newton(fname,x);

xx = linspace(-2, 4)';
yy = linspace(-2,  4);

for k = 1:100
    aux = xx(k);
    for j =1:100
        Z(k,j)= feval(fname, [aux, yy(j)]');
    end
end

[XX,YY] = meshgrid(xx,yy);

x1 = W(1,:); x2 = W(2,:);

contour(XX,YY,Z, 1000)
hold on

plot(x1,x2,'or', x1, x2,'--b','Linewidth',3)
title(' Convergencia del método','Fontsize',16)
