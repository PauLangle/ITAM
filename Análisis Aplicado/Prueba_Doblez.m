%Script file: Prueba_Doblez
% Prueba la región de confianza y la reducción del modelo.
%
% Análisis Aplicado
%   ITAM
%  21 de septiembre de 2020.
%------------------------------------------------------

close all

fname = 'rosenbrock';   % nombre de la función
%fname = 'Shubert';
xs = [1.5 2.5]';           % punto anclado
f = feval(fname,xs);        % valor de la función en el punto
g = gradiente(fname,xs);    % gradiente en el punto
B = hessian(fname,xs);      % matriz hessiana en el punto
ps = - B\g;                 % Dirección de Newton en el punto
%--------------------------------------------------------------------
% dirección de Newton y punto de Cauchy
t = linspace(0,1,30)';
pc = -((g'*g)/(g'*B*g))*g;
d = ps - pc;
pp = zeros(2,30); qq = zeros(2,30);
for k = 1:30
    pp(:,k) = xs + t(k)*ps;
    qq(:,k) = xs +t(k)*pc;
end
plot(pp(1,:),pp(2,:),'g', qq(1,:),qq(2,:),'m','linewidth',3)
legend('Dirección de Newton', 'Punto de Cauchy')
hold on
%-------------------------------------------------------

iter =  0;                  % iteraciones internas

P = [];                     % Se guardan las soluciones

Delta = 0.5;
% Parte iterativa
while(iter < 5)
    [pn] = doblez(B,g,Delta);   % Doblez
     P = [P pn];   % de guarda el punto  de la región de confianza
     Delta = Delta/2;
    iter = iter + 1;                    % se incrementan las iteraciones
    theta = linspace(0,2*pi,50)';
    for i = 1:50
        cir(1:2, i) = Delta*[cos(theta(i)), sin(theta(i))]';
    end
    plot( cir(1,:)+xs(1), cir(2,:)+xs(2),'k','Linewidth',2)
    axis equal
    hold on
end

%---------------------------------------------------------------------------------------------
% Graficación

plot(P(1,:)+xs(1),P(2,:)+xs(2),'dr', P(1,:)+xs(1), P(2,:)+xs(2), 'b','LineWidth',3)
title(' Aproximación a Región de Confianza/Doblez/Rosenbrock', 'FontSize',16)
hold on
plot(xs(1),xs(2), 'dk', 'Linewidth',3)

