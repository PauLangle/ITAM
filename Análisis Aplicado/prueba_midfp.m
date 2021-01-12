%   Prueba_midfp.m
%----------------------------------------------------------------
% Equipo:
% Pablo Martines CU: 165680
% Paula Langle CU:162574
% Ricardo Zuñiga CU: 157116
%
% fecha: 2 de diciembre 2020
%-----------------------------------------------------------------
% Este script resuelve el problema de Mínimos Cuadrados no Lineales
% para la función logsitica c(t)=x(3)/(1+x(2)*exp(-x(1)*t)) con los datos
% de cosechas de trigo por tonelada y hectarea en los ultimos 24 anios. 
% El objetivo es encontrar los parametros optimos (x1 x2 x3) para el problema de la
% cosecha utilizando el metodo DFP a la inversa con busqueda de linea.
% ------------------------------------------------------------------------
% Primero desplegamos la tabla con la solución

close all
f='ftrigo';
x0=[0,1,30]';

[x,Z,iter,sol]=dfp(f,x0);

a=x(1);
b=x(2);
c=x(3);
Norma=Z(end);
fprintf("\n----------------------------------------------------------\n");
resultado=table(a,b,c,Norma,iter);
disp(resultado);
%----------------------------------------------------------------------
% Ahora preparamos el codigo para la graficacion de la solución
% Primero pasamos los datos

t=[1:24]'; %Línea del tiempo
%Datos de la cosecha por año, por tonelada por hectarea
y=[11.72 13.38 14.10 13.87 14.80 15.58 14.36 16.30 16.91 18.16];
y=[y 18.43 18.70 20.46 19.16 20.01 22.41 21.21 22.81 23.97];
y=[y 23.27 23.80 25.59 24.93 26.59];
y=y'; %Vector columna

plot(t,y,'rd','Linewidth',3)
title('Función solución al problema de Mínimos Cuadrados no Lineales','Fontsize',12)
hold on

%xin=[0.1 1 30]';
c1=x(3);

tt=linspace(1,24,200)'; %Partición

% Gráfica con mis valores iniciales (color verde)
% función c(t)=x(3)/(1+x(2)*exp(-x(1)*t));
vv=c1./(1+b.*exp(-a*tt));
plot(tt,vv,'b','Linewidth',3);
xlabel('Años') 
ylabel('Cosecha (tonelada x hectarea)') 
legend({'Datos reales','función solución'},'Location','northwest')
hold off

%-------------------------------------------------------------------------------
% Gráfica en tres dimensiones para la solución

a=[sol(1)];
b=[sol(2)];
c=[sol(3)];

for i=1:3:size(sol)-5
    a=[a;sol(i+3)];
end

for i=2:3:size(sol)-4
    b=[b;sol(i+3)];
end

for i=3:3:size(sol)-3
    c=[c;sol(i+3)];
end


%plot3(a,b,c)
