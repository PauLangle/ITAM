% Script: calor2.m
% Distribuci�n de temperatura en una placa met�lica
% Se resuelve un sistema lineal Ax=b
% Con A sim�trica positiva definida y usamos gradiente conjugado
%
% An�lisis Aplicado
% 5 de octubre
% --------------------------------------------------------------

m=30;
[A,b]=matcalor(m);
[x,k]=migc(A,b);

u=(1:m)';
v=(1:m)';
Z=zeros(m); 

% x es un vector de m*m
% x(1:m) temperatura en el primer rengl�n de mi maya

for k=1:m
    Z(k,1:m)=x((k-1)*m+1:m*k)'; %En Z separo los valores de la x para poder 
                               % armar la maya con los valores ya resueltos
                               % de x. 
end

%La graficaci�n
[U,V]=meshgrid(u,v);
surf(U,V,Z)
title('Distribuci�n de Calor en una Placa Met�lica','Fontsize',16)