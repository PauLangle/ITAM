function [fx] = fcuad2(x)
% Funcion cuadratica f: R^n --> R

% In
% x .- vector de longitud n
%Out
%fx .- número real

a1 = 100; a2 = 10; %Estas son las lambdas
A = [a1 0 ; 0 a2]';

y = [x(1)-1; x(2)-2];

 fx = (1/2)*y'*A*y +100;
 
end