%Script: pruebamigc.m
m=30;
[A,b]=matcalor(m);
%A es nXn con n=m*m
%b es un vector nX1
tic;
[x,iter]=migc(A,b);
toc

fprintf('Iteraciones: %3.0f\n',iter);

tic;
[x1,flag,relres,iter1]=pcg(A,b,1e-08); %Método de Matlab
toc

fprintf('Iteraciones método 2: %3.0f\n',iter1);

%norm(x-x1)  %Vemos la diferencia entre los vectores resultado 

% migc.m se usa en las siguientes funciones:
% doblez.m -> pN=-B\g cambia a pN=migc(B,-g)
% miregion.m

% Rosenbrock100.m 
% miregion.m en Rosenbrock100.m con x0=3.5*ones(100,1)
% Resultados:
% iter= 19 o 22
% x* =ones(100,1)
