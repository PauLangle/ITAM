function [x,iter] = descensonewton(fname,x)
% Método de búsqueda  usando dirección de Newton.

% Análisis Aplicado
% ITAM
% 26 de agosto de 2020.
% Dr. Zeferino Parada

%In
% fname.- cadena de caracteres con el nombre de la función a minimizar.
% x.- vector n-dimensional.
% Out
% x.- vector n-dimensional con la aproximación al mínimo local.
% iter.- contador con el número final de iteraciones externas.
%
% Ll�mala as�: [x,iter]=descensonewton('rosenbrock',[6,8]')
%
% Se requiere aproximaciones al gradiente en la función
%  gradiente.m con el llamado:
% gfx = gradiente (fname,x);
% La matriz hessiana se aproxima con la función:
% H = hessian(fname,x);

% parámetros
tol     = 1.e-05; % tolerancia para la norma del gradiente.
tolp = 1.e-08;    % tolerancia para tamaño de paso en cada iteración
maxiter = 100;    % número máximo de iteraciones externas permitidas
maxjter = 20;
c1 = 0.1;
% valores iniciales
iter = 0;        % contador para las iteraciones externas

g     = gradiente(fname,x);
ng    = norm(g);
paso = 1;
disp('iter         || Df(x_k) ||'  )
disp('-----------------------------------')
%-----------------------------------------------------------------------
% parte iterativa
while ( ng > tol && iter < maxiter && paso>tol )
     H = hessian(fname,x);    % dirección de Newton
     p = - H \ g;
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
    paso = alfa*norm(p);
     x = x + alfa*p;
      g = gradiente(fname,x);
      ng = norm(g);
      iter = iter + 1; 
      disp(sprintf('%2.0f %2.8f  %2.8f',iter, ng,  paso))
end

end
