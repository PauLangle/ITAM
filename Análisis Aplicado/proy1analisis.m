%----------------------------------------------------------------
%Equipo:
% Pablo Martines CU: 165680
% Paula Langle CU:162574
% Ricardo Zuñiga CU: 157116
%-----------------------------------------------------------------

%Este script prueba el método del doblez (con la aceptación del paso) para aproximar el mínimo local de
%la función de Rosenbrock.

%Definimos nuestros parámetros de entrada
f ='rosenbrock';
x0 = [3.5,4.5];

%Medimos el tiempo de ejecución del método de doblez
tic;
[xf, iter,ng] = miregion(f, x0);
toc;

%Desplegamos la tabla con los resultados:
fprintf('Resultado:\n')
x1=xf(1);
x2=xf(2);
norma=ng;
fprintf('\n')
resultado=table(iter,x1,x2,norma);
disp(resultado);
