% Script: Prueba_maxdesc.m
close all

fname = 'fcuad2';

xs = [1 2]';
%Graficamos  la función cuadrática alrededor de su mínimo local
n = 75;
x = linspace(-1,3,n)'; y = linspace(0,4,n);


for k = 1:n
    aa = x(k);
    for j = 1:n
        v = [aa,y(j)]';
        if (norm(v-xs) <= 2)
             fv = feval(fname,[aa y(j)]');
             plot3(aa,y(j),fv,'ob',aa,y(j),0,'xy')
             hold on
        end
    end
end
plot3(xs(1),xs(2),0,'dr', xs(1),xs(2),100+5, 'dr','Linewidth',3)
xlabel('Eje X1', 'Fontsize', 16)
ylabel('Eje X2','Fontsize',16)
title(' f(x) = 100*(x_1 -1)^2 + 10*(x_2 -2)^2 + 100','Fontsize',16)
hold on

w = linspace(0,100,200);
for k = 1:200
   plot3(xs(1),xs(2),w(k) +5,'xr')
   hold on
end

x = [2 8]';
[x,iter,W] = Maxdesc(fname,x);

