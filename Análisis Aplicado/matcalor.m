function [A,b] = matcalor(n);
% Esta función genera la matriz y el lado derecho
% del problema de calor de una pieza rectangular.
% en una malla rectangular de nxn.
% La matriz A es simétrica positiva definida de orden (n*n) y es la
% matriz de Poisson.

% el calor en la parte superior es de 100 grados,
% el calor en la parte inferior es de 25 grados,
% el calor en las partes laterales es de 50 grados,
% el número de puntos donde se debe aproximar el calor
% es n^2.
% Calor en el punto k-ésimo es el promedio del calor
% de sus cuatro vecinos, al norte, sur, este y oeste.
% La numeración de los nodos va por renglones
% de izquierda a derecha iniciando en la esquina 
% superior izquierda.
% Ejemplo:
%         1       2       3      4
%         x ------x ----- x -----x
%         |       |       |      |
%         |       |       |      |
%         x ------x-------x------x
%         5       6       7      8

% In
%  n.- número natural mayor a 1.
% Out
%  A.- matriz n2xn2 con n2 = n*n.
%      El renglón i-ésimo de A es la ecuación del
%      punto i-ésimo en el arreglo rectangular.
%  b. - vector columna de orden n2.
%       la entrada i-ésima es el lado derecho de
%       la ecuación i-ésima.
%
%----------------------------------------------------
%
%
%  Análisis Aplicado
% ITAM
%  30 de septiembre de 2020.
%---------------------------------------------------

% inicio
n2 = n^2;
A  = 4*eye(n2);

% subdiagonal
for k = 2:n2
    A(k,k-1)= -1.0;
end
% depurando la subdiagonal
for k = 1:n-1
    A(k*n +1, k*n) = 0.0;
end

% subdiagonal -n
for k = 1:n2-n
    A(n+k,k) = -1.0;
end

A = A + (tril(A,-1))';

% lado derecho
b = zeros(n2,1);

b(1) = 150.0;
b(n) = 150.0;

for k = 2:n-1
   b(k) = 100.0;
end

for k = (n+1):n:(n*(n-2) +1)
   b(k) = 50.0;
end

for k = 2*n:n:(n*(n-2))
   b(k) = 50.0;
end

b(n2-n +1 ) = 75.0;
b(n2) = 75.0;

for k = (n2 -n +2): (n2-1)
   b(k) = 25.0;
end
end