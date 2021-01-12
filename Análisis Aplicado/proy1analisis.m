%----------------------------------------------------------------
%Equipo:
% Pablo Martines CU: 165680
% Paula Langle CU:162574
% Ricardo Zu�iga CU: 157116
%-----------------------------------------------------------------

%Este script prueba el m�todo del doblez (con la aceptaci�n del paso) para aproximar el m�nimo local de
%la funci�n de Rosenbrock.

%Definimos nuestros par�metros de entrada
f ='rosenbrock';
x0 = [3.5,4.5];

%Medimos el tiempo de ejecuci�n del m�todo de doblez
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
