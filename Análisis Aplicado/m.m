function fx = m(p,B,g,x)
%Función correspondiente al modelo cuadrático
%In
% p - vector de longitud 2
% B - matriz Hessiana evaluada en un punto
% g - gradiente de la función evaluada en un punto
%Out
% fx.- número real con el valor de la función.

fx=(1/2)*p'*B*p+g'*p+feval('rosenbrock',x);

