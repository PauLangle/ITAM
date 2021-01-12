function [x,iter] = MetodoBL(fname,x)
% M�todo de b�squeda de l�nea con la primer condici�n de Wolfe
% usando m�ximo descenso.
% 12 de agost de 2020.
% ITAM
% An�lisis Aplicado
% Dr. Zeferino Parada
%
%In
% fname.- cadena de caracteres con el nombre de la funci�n a minimizar.
% x.- vector n-dimensional.
% Out
% x.- vector n-dimensional con la aproximaci�n al m�nimo local.
% iter.- contador con el n�mero final de iteraciones externas.
%
% Se requiere aproximaciones al gradiente en la funci�n
%  gradiente.m con el llamado:
% g = gradiente (fname,x);
%------------------------------------------------------------------------

% par�metros
tol     = 1.e-04; % tolerancia para la norma del gradiente.
c1      = 0.1;    % valor para la condici�n de Wolfe
maxiter = 50;    % n�mero m�ximo de iteraciones externas permitidas
maxjter = 30;     % n�mero m�ximo de iteraciones internas de BL permitidas

% valores iniciales
iter = 0;        % contador para las iteraciones externas
fx = feval(fname,x);
g     = gradiente(fname,x);
ng    = norm(g);
p = 1; s = -1;
% parte iterativa
while ( ng > tol && iter < maxiter &&  norm(p)> 1.e-04)
     %p = -g;             % m�ximo descenso
     H = hessian(fname,x);
     p = -H\g;  
     
    %-----------------------------------------
    % B�squeda de l�nea
    alfa  = 1;              % paso completo
    xt = x + alfa*p;        % primer punto de prueba
    f  = feval(fname,x);  % valor de la funci�n
    f1 = feval(fname, xt);% valor de la funci�n en el punto de prueba
    s  = p'*g;            % derivada direccional
    jter = 0;             % iteraciones internas
    
    if ( abs(s) < 1.e-06)
        disp('No existe suficiente descenso  ')
        disp(sprintf('%2.0f  %2.10f',iter, s   ))
        iter = maxiter; 
    end
    while( (f1>f+alfa*c1*s) && jter < maxjter)  % b�squeda de l�nea
       alfa = alfa/2;
       xt = x + alfa*p;
       f1 = feval(fname, xt);% valor de la funci�n en el punto de prueba
       jter = jter +1;
    end   
    % Fin de b�squeda de l�nea
    %----------------------------------------------------------------------
    
     %--------------------------------------------------------------
     % Graficaci�n
       t = linspace(0,1,50)';
       ft = zeros(50,1); rt = zeros(50,1);
       for k = 1:50        
          ft(k) = feval(fname, x+t(k)*p);  % funci�n 
          rt(k) = fx + t(k)*(c1*p'*g);
       end
       fx = feval(fname,x+alfa*p);
       plot(t,ft,'--b',t,rt,'--m',alfa,fx,'dr', 'LineWidth',3)
       title('Gr�fica de b�squeda de l�nea','Fontsize',16)
       xlabel('EJE  T','Fontsize',16)
       ylabel(' f(x + tp)','Fontsize',16  )
       legend('f(x)','recta','punto')
        pause
      % Fin de graficaci�n
      %-----------------------------------------------------  
      % Actualizaci�n de valores
     x = x + alfa*p;  
     fx = feval(fname,x); 
     g = gradiente(fname,x);
     iter = iter + 1;
end       % Fin del while principal


end


