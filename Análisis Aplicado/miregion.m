function [xf,iter,ng] = miregion(fname,x0)
% Metodo de region de confianza para aproximar un minimo
% de f : Rn -> R dos veces continuamente diferenciable
% In
% f cadena de caracteres con el codigo en Matlab de la funcion a minimizar.
% x0 vector columna de dimension n con el punto inicial.
% Out
% xf vector columna de dimension n con la aproximacion final.
% iter numero de iteraciones que se realizaron.
%----------------------------------------------------------------
%Equipo:
% Pablo Martines CU: 165680
% Paula Langle CU:162574
% Ricardo Zuñiga CU: 157116
%-----------------------------------------------------------------
% Criterios para parar
%||grad(f(xk))||<=tol o iter >=maxiter o jregion>=maxregion
% Parametro maxregion es el numero maximo que permanece
% region de confianza en un solo punto.
% Parametros:
deltamin = 1.e-04;
deltamax = 5;
nu = 0.25;
maxiter = 100;
maxregion = 10;
tol = 1.e-06;
ros=0;

xs =x0';           %Transponemos nuestro punto inicial
delta=1;           %Ponemos una delta inicial
iter=0;            %Inicializamos el número de iteraciones
ng=0;              %Inicializamos el valor de la norma

xs2=xs;            %Variable auxiliar que nos ayudará a guardar el paso anterior
jregion=0;         %Inicializamos el contador que nos indica cuántas veces rechazamos un paso. 

g = gradiente(fname,xs);    % Gradiente en el punto inicial.
ng=norm(g);

while(ng>tol && iter<maxiter && jregion<maxregion)
    iter=iter+1;
    B = hessian(fname,xs);  % Matriz hessiana en el punto
    ps=doblez(B,g,delta);   %Calculamos el punto que minimiza el modelo cuadrático
    
    %Ahora haremos la aceptación del paso:
    ros=(feval(fname,xs)-feval(fname,xs+ps))/(-(1/2)*ps'*B*ps-g'*ps);
    
    %Guardamos el paso anterior
     xs2=xs;
    
    if ros>=(1-nu)
        %Aceptamos el paso, pero los podemos dar más grandes.
        %Checamos que la delta no pase de DeltaMax
        if delta*2<=deltamax
            delta=2*delta;
        else
            delta=deltamax;
            %jregion=jregion+1;
        end
        
        xs=xs+ps;
    else
        if nu<=ros && ros<(1-nu)
            %Se acepta el paso y no es necesario mover delta
            xs=xs+ps;
        else
            %No se acepta el paso, hay que hacerlo más chico.
            %Cuida que no se pase de la cota inferior.
            
            if 0.5*delta>=deltamin
                delta=0.5*delta;
            else
                delta=deltamin;
                jregion=jregion+1; %Sumamos al contador que nos indica cuántas veces
            end                    % rechazamos el paso k. 
        end
    end
    
    
    % Si ya nos movimos tenemos que reiniciar el contador de jregion para
    % el nuevo punto, Xs2 es el paso anterior. 
    if xs~=xs2
        jregion=0;
    end
    
     %Volvemos a calcular el gradiente
     g = gradiente(fname,xs);  
     ng=norm(g);   
     
     
end  
 

%Preparamos la salida del método

if jregion>=maxregion
    fprintf('Salimos del método porque no nos estamos moviendo del punto.\n')
    xf=xs;
end
 

if iter>=maxiter
        fprintf('Se sobrepasó el número de iteraciones permitidas.\n')
end


if norm(g)<=tol
        fprintf('Encontramos una solución.')
        xf=xs;
end
        
   


   
