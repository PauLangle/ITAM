% Script: calor2.m
% Distribución de temperatura en una placa metálica
% Se resuelve un sistema lineal Ax=b
% Con A simétrica positiva definida y usamos gradiente conjugado
%
% Análisis Aplicado
% 5 de octubre
% --------------------------------------------------------------

m=30;
[A,b]=matcalor(m);
[x,k]=migc(A,b);

u=(1:m)';
v=(1:m)';
Z=zeros(m); 

% x es un vector de m*m
% x(1:m) temperatura en el primer renglón de mi maya

for k=1:m
    Z(k,1:m)=x((k-1)*m+1:m*k)'; %En Z separo los valores de la x para poder 
                               % armar la maya con los valores ya resueltos
                               % de x. 
end

%La graficación
[U,V]=meshgrid(u,v);
surf(U,V,Z)
title('Distribución de Calor en una Placa Metálica','Fontsize',16)