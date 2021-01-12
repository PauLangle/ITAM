function [x,iter,W] = maxdesc(fname,x)
% Método para aproximar mínimos locales para
% fname: R^n --> R  dos veces continuamente diferenciable,
% usando el negativo del gradiente.

% Análisis Aplicado
% ITAM
% 24 de agosto de 2020.
% Dr. Zeferino Parada

%In
% fname.- cadena de caracteres con el nombre de la función a minimizar.
% x.- vector n-dimensional.
% Out
% x.- vector n-dimensional con la aproximación al mínimo local.
% iter.- contador con el número final de iteraciones externas.
%
% Se requiere aproximaciones al gradiente en la función
%  gradiente.m con el llamado:
% gfx = gradiente (fname,x);

%---------------------------------------------------------------------
% parámetros
tol     = 1.e-05; % tolerancia para la norma del gradiente.
maxiter = 100;    % número máximo de iteraciones externas permitidas
maxjter = 50;
c1 = 0.1;
% valores iniciales
iter = 0;        % contador para las iteraciones externas

n = length(x);
g     = gradiente(fname,x);
ng    = norm(g);
W = [x];  % matriz de graficación
%-----------------------------------------------------------------------
% parte iterativa
while ( ng > tol && iter < maxiter)
    
    p = -g;  % máximo descenso
     
    %---------------------------------------------
    alfa  = 1;              % paso completo
    xt = x + alfa*p;
    % busqueda de línea
    f  = feval(fname,x);  % valor de la función
    f1 = feval(fname, xt);% valor de la función en el punto de prueba
    s  = p'*g;            % derivada direccional
    jter = 0;             % iteraciones internas
    while( (f1>f+alfa*c1*s) && jter < maxjter)  % búsqueda de línea
       alfa = alfa/2;
       xt = x + alfa*p;
       f1 = feval(fname, xt);% valor de la función en el punto de prueba
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



