function [x,iter] = MetodoBL(fname,x)
% Método de búsqueda de línea con la primer condición de Wolfe
% usando máximo descenso.
% 12 de agost de 2020.
% ITAM
% Análisis Aplicado
% Dr. Zeferino Parada
%
%In
% fname.- cadena de caracteres con el nombre de la función a minimizar.
% x.- vector n-dimensional.
% Out
% x.- vector n-dimensional con la aproximación al mínimo local.
% iter.- contador con el número final de iteraciones externas.
%
% Se requiere aproximaciones al gradiente en la función
%  gradiente.m con el llamado:
% g = gradiente (fname,x);
%------------------------------------------------------------------------

% parámetros
tol     = 1.e-04; % tolerancia para la norma del gradiente.
c1      = 0.1;    % valor para la condición de Wolfe
maxiter = 50;    % número máximo de iteraciones externas permitidas
maxjter = 30;     % número máximo de iteraciones internas de BL permitidas

% valores iniciales
iter = 0;        % contador para las iteraciones externas
fx = feval(fname,x);
g     = gradiente(fname,x);
ng    = norm(g);
p = 1; s = -1;
% parte iterativa
while ( ng > tol && iter < maxiter &&  norm(p)> 1.e-04)
     %p = -g;             % máximo descenso
     H = hessian(fname,x);
     p = -H\g;  
     
    %-----------------------------------------
    % Búsqueda de línea
    alfa  = 1;              % paso completo
    xt = x + alfa*p;        % primer punto de prueba
    f  = feval(fname,x);  % valor de la función
    f1 = feval(fname, xt);% valor de la función en el punto de prueba
    s  = p'*g;            % derivada direccional
    jter = 0;             % iteraciones internas
    
    if ( abs(s) < 1.e-06)
        disp('No existe suficiente descenso  ')
        disp(sprintf('%2.0f  %2.10f',iter, s   ))
        iter = maxiter; 
    end
    while( (f1>f+alfa*c1*s) && jter < maxjter)  % búsqueda de línea
       alfa = alfa/2;
       xt = x + alfa*p;
       f1 = feval(fname, xt);% valor de la función en el punto de prueba
       jter = jter +1;
    end   
    % Fin de búsqueda de línea
    %----------------------------------------------------------------------
    
     %--------------------------------------------------------------
     % Graficación
       t = linspace(0,1,50)';
       ft = zeros(50,1); rt = zeros(50,1);
       for k = 1:50        
          ft(k) = feval(fname, x+t(k)*p);  % función 
          rt(k) = fx + t(k)*(c1*p'*g);
       end
       fx = feval(fname,x+alfa*p);
       plot(t,ft,'--b',t,rt,'--m',alfa,fx,'dr', 'LineWidth',3)
       title('Gráfica de búsqueda de línea','Fontsize',16)
       xlabel('EJE  T','Fontsize',16)
       ylabel(' f(x + tp)','Fontsize',16  )
       legend('f(x)','recta','punto')
        pause
      % Fin de graficación
      %-----------------------------------------------------  
      % Actualización de valores
     x = x + alfa*p;  
     fx = feval(fname,x); 
     g = gradiente(fname,x);
     iter = iter + 1;
end       % Fin del while principal


end


