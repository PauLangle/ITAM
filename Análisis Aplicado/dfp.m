function [x,Z,iter,sol] = dfp(fname,x0)
% Optimizacion con busqueda de linea
% utilizando el metodo DFP a la inversa.
%
% Parametros de entrada
% fname : nombre de la funcion a optimizar.
% x0 :  punto inicial del proceso de minimización
%
% Parametros de salida
%
% x : aproximacion a un minimo local de fname.
% Z: vector con las normas del gradiente en cada paso
% iter : iteraciones necesarias para llegar al optimo
% sol: vector que contiene la solucion en la iteracion K-esima

%----------------------------------------------------------------------
tol = 1.e-05;                    % tolerancia de la norma del gradiente
x   = x0;                        % punto inicial
n   = length(x);       
B   = eye(n);                    % matriz inicial para el metodo
gx  = gradiente(fname,x);        % gradiente en el punto inicial

iter    = 0;                     % contador de iteraciones
maxiter = 100;                   % máximo número de iteraciones          
maxjter = 10;                    % número máximo de iteraciones en busqueda de linea
Z = [norm(gx)];                  % vector que guarda las normas del gradiente
sol=[x];                         % vector que va guardando la solucion en la iteracion kesima
paso=1.0;

c1 =1e-02;                      % Valor de la constante para condicion de Wolf

while ( norm(gx) > tol  && iter < maxiter && paso > 1.e-05)
    % Calculo de direccion de descenso sin resolver el sistema lineal
    p =- B*gx;  
    
    % Comenzamos con busqueda de lineal
    j = 0;             %El contador jter y alfa se deben de reiniciar en cada iteracion
    alfa = 1.0;           
    pend = c1*(gx'*p);      %Calculo de la pendiente en condicion de wolf
    fx = feval(fname,x);  %Evaluamos la funcion en el punto x
    xn = x +alfa*p;
    fxn=feval(fname,xn); %Evaluamos la funcion dando el nuevo paso
    
    while ( (fxn >= fx + alfa*pend ) && j < maxjter ) %Checamos que se cumpla Wolf
        alfa = (1/2)*alfa;  %Si es necesario, achicamos el paso
        xn = x +alfa*p;     % Calculamos el nuevo punto
        fxn=feval(fname,xn); %Actualizamos el valor de f en el nuevo paso
        j = j + 1;
    end
            
    xp = x + alfa*p; %Una vez aceptado el paso, avanzamos
    
    gxp = gradiente(fname,xp); %Calculamos el gradiente de la funcion en el nuevo punto
    
    % Nos preparamos para actualizar la matriz B
    s = xp - x ;
    y = gxp - gx;
    gama = 1/(s'*y);
    
    if (s'*y > 0)  % Checamos que se cumpla la condicion necesaria para relizar la actualizacion
        flag = 1;  % Bandera para saber si actualizamos con el metodo deseado
        B = (eye(n) - gama*s*y')*B*(eye(n)-gama*y*s') + gama*(s*s'); %Actualizamos la matriz
    else
        flag = 0;  %Bandera para saber si no se pudo actualizar con el metodo deseado
        B = eye(n); %En dicho caso, volvemos a la matriz identidad
    end

    % Actualizamos el valor de nuestras variables
    paso = norm(xp - x);
    x = xp;
    gx = gxp;
    Z = [Z;norm(gx)];
    sol=[sol;x];
    iter = iter+1;
    disp(sprintf('%2.0f   %2.14f   %1.0f', iter,  norm(gx), flag ))
end