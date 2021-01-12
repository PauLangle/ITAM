%Script cosecha.m

%Grafica la gr�fica de la cosecha

t=[1:24]'; %L�nea del tiempo
%Datos de la cosecha por a�o, por tonelada por hectarea
y=[11.72 13.38 14.10 13.87 14.80 15.58 14.36 16.30 16.91 18.16];
y=[y 18.43 18.70 20.46 19.16 20.01 22.41 21.21 22.81 23.97];
y=[y 23.27 23.80 25.59 24.93 26.59];
y=y'; %Vector columna

plot(t,y,'rd','Linewidth',3)
title('Valores de la cosecha de trigo','Fontsize',16)
hold on

xin=[0.1 1 30]';
a=0.1; b=1; c1=30;

tt=linspace(1,24,200)'; %Partici�n

% Gr�fica con mis valores iniciales (color verde)
% funci�n c(t)=x(3)/(1+x(2)*exp(-x(1)*t));
vv=c1./(1+b.*exp(-a*tt));
plot(tt,vv,'g','Linewidth',3);

%Usando la funci�n de Matlab, gr�fica �ptima (color azul)
[x,fx,exitflag,output]=fminunc('ftrigo',xin);

a=x(1); b=x(2); c1=x(3);

vv=c1./(1+b.*exp(-a*tt));

plot(tt,vv,'b','Linewidth',3)
