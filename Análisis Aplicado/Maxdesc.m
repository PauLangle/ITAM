function [x,iter,W] = maxdesc(fname,x)
% M�todo para aproximar m�nimos locales para
% fname: R^n --> R  dos veces continuamente diferenciable,
% usando el negativo del gradiente.

% An�lisis Aplicado
% ITAM
% 24 de agosto de 2020.
% Dr. Zeferino Parada

%In
% fname.- cadena de caracteres con el nombre de la funci�n a minimizar.
% x.- vector n-dimensional.
% Out
% x.- vector n-dimensional con la aproximaci�n al m�nimo local.
% iter.- contador con el n�mero final de iteraciones externas.
%
% Se requiere aproximaciones al gradiente en la funci�n
%  gradiente.m con el llamado:
% gfx = gradiente (fname,x);

%---------------------------------------------------------------------
% par�metros
tol     = 1.e-05; % tolerancia para la norma del gradiente.
maxiter = 100;    % n�mero m�ximo de iteraciones externas permitidas
maxjter = 50;
c1 = 0.1;
% valores iniciales
iter = 0;        % contador para las iteraciones externas

n = length(x);
g     = gradiente(fname,x);
ng    = norm(g);
W = [x];  % matriz de graficaci�n
%-----------------------------------------------------------------------
% parte iterativa
while ( ng > tol && iter < maxiter)
    
    p = -g;  % m�ximo descenso
     
    %---------------------------------------------
    alfa  = 1;              % paso completo
    xt = x + alfa*p;
    % busqueda de l�nea
    f  = feval(fname,x);  % valor de la funci�n
    f1 = feval(fname, xt);% valor de la funci�n en el punto de prueba
    s  = p'*g;            % derivada direccional
    jter = 0;             % iteraciones internas
    while( (f1>f+alfa*c1*s) && jter < maxjter)  % b�squeda de l�nea
       alfa = alfa/2;
       xt = x + alfa*p;
       f1 = feval(fname, xt);% valor de la funci�n en el punto de prueba
       jter = jter +1;
    end   
    %--------------------------------------
     x = x + alfa*p;
     
     W = [W x];
      g = gradiente(fname,x);
      ng = norm(g);
      iter = iter + 1; 
      disp(sprintf('%2.0f % 2.10f',iter, ng  )    )
end



