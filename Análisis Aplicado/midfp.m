function [x, k, Z]= midfp(f,x0)
% Metodo de busqueda de linea para f: R^n --> R con actualizacion DFP
%
% Proyecto Analisis Aplicado
% ITAM
% 24 de noviembre 2020
%

% parametros
tol = 1e-05;  % tolerancia  de la norma gx
kmax = 100;   % numero maximo de iteraciones
c1 = 1.e-04;   % primer condicion de Wolfe
k = 0;        % contador de iterciones
Z = [];
% Valores iniciales
 n = length(x0);
 B = eye(n);         % matriz inicial
 x = x0;
 gx = gradiente(f,x);  % gradiente de f en x
 %B = norm(gx)*eye(n);
 paso = 1.0;
 
 while( norm(gx)> tol &&  k < kmax && paso > 1.e-05)
    % direccion de descenso: solucion de  B*p = -gx 
    p = -B\gx;
    
    % Busqueda de linea/ primera condicion de Wolfe
    fx = feval(f,x);
    jmax = 10; j = 0;
    pen = c1*(gx'*p);
    alfa = 1;
    xn = x +alfa*p;
    fxn = feval(f,xn);
    
    %Checamos la aceptacion del paso con Wolf
    while( (fxn > fx + alfa*pen)  && j < jmax  )
       alfa = alfa/2;
       xn = x +alfa*p; 
       fxn = feval(f,xn);
       j = j+1;
    end
    
    
    xp = x + alfa*p;
    gxp = gradiente(f,xp);
    
    % Nos preparamos para actualizar la matriz
    s = xp-x;
    y = gxp-gx;
    
    %Actualizacion de DFP en caso de que NO tengamos que resolver un sistema de ecuaciones
     if (y'*s > 0)
        %Si se generó una dirección de descenso
        flag = 1;
        gama = 1 / (s' * y);
        %B = (eye(n) - gamma * y * s') * B * (eye(n) - gamma * s * y') + gamma * (y * y');
        B = (eye(n) - gama*s*y')*B*(eye(n)-gama*y*s') + gama*s*s';
    else
        %Si no se generó una dirección de descenso
        flag = 0;
        B = eye(n);
    end

% Actualizacion de DFP en caso de que SI tengamos que resolver un sistema de ecuaciones
%     gama = 1 / (s' * y);
%     B = (eye(n) - gamma * y * s') * B * (eye(n) - gamma * s * y') + gamma * (y * y');
%     flag = 1;
%     
%     if (cond(B) > 1e+04)
%         B = eye(n);
%         flag = 0;
%     end

        
    
    % Actualizar valores
     paso = norm(xp - x);
     x = xp;
     gx = gxp;
     k = k+1;
     Z = [Z; norm(gx)];
     fprintf('%2.0f %2.8f %2.0f %2.8f\n', k, norm(gx), flag, alfa) 
 end

end