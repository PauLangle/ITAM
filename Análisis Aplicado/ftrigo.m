function [fx]=ftrigo(x)
%Función f: R^n --> R con la cosecha del trigo

t=[1:24]'; %Línea del tiempo
%Datos de la cosecha por año, por tonelada por hectarea
y=[11.72 13.38 14.10 13.87 14.80 15.58 14.36 16.30 16.91 18.16];
y=[y 18.43 18.70 20.46 19.16 20.01 22.41 21.21 22.81 23.97];
y=[y 23.27 23.80 25.59 24.93 26.59];
y=y'; %Vector columna

% función c(t)=x(3)/(1+x(2)*exp(-x(1)*t));

fx=0;
a=x(1); b=x(2); c1=x(3); %Por eficiencia sacamos los valores
for k=1:24
    nv= c1/(1+b*exp(-a*t(k))); %Evaluas la funcion
    nv=nv-y(k);  %Sacas el residuo
    nv=nv*nv;    % Elevas al cuadrado
    fx=fx+nv;
end

fx=fx/2;
end