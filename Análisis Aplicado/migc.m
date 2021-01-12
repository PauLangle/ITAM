function [x,iter]=migc(A,b)
%M�todo del gradiente conjugado para resolver el sistema lineal 
% Ax=b 
%donde A es una matriz nxn sim�trica positiva definida
% versi�n #1
%
% An�lisis Aplicado 
% ITAM 
% 5 de octubre 2020
%---------------------------------------------------------------

n=length(b);    %dimensi�n del problema
x=zeros(n,1);   %punto inicial

r=A*x-b;        %primer residual
p=-r;           %primera direcci�n
itermax=10*n;
iter=0;         %contador de iteraciones
tolr=1e-08;     % tolerancia a la norma del residual
            

while(norm(r)>=tolr && iter<itermax)
    residual=r'*r;
    alfa=(residual)/(p'*A*p);  %Calcular alfa
    x=x+alfa*p;                 %Nuevo punto
    r=A*x-b;                    %Nuevo residual
    beta=(r'*r)/(residual);     %Calcular beta
    p=-r+beta*p;                %Nueva direcci�n
    iter=iter+1;                %Actualizamos el contador
end

end
